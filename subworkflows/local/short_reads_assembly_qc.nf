include { DECONTAMINATE_ASSEMBLIES } from '../../subworkflows/ebi-metagenomics/decontaminate_assemblies/main'
include { SEQKIT_SEQ               } from '../../modules/nf-core/seqkit/seq/main'

workflow SHORT_READS_ASSEMBLY_QC {

    take:
    assembly               // [ val(meta), path(assembly_fasta) ]

    main:
    ch_versions = channel.empty()

    /* Len filter using the parameter "short_reads_min_contig_length" */
    SEQKIT_SEQ(assembly)
    ch_versions = ch_versions.mix(SEQKIT_SEQ.out.versions)

    /***************************************************************************/
    /* Remove PhiX, human and host contigs (each step runs only if requested)  */
    /***************************************************************************/
    DECONTAMINATE_ASSEMBLIES(SEQKIT_SEQ.out.fastx)
    ch_versions = ch_versions.mix(DECONTAMINATE_ASSEMBLIES.out.versions)

    def cleaned_contigs = DECONTAMINATE_ASSEMBLIES.out.cleaned_contigs

    def contaminated_contigs_by_step = DECONTAMINATE_ASSEMBLIES.out.phix_contaminated_contigs_tsv.map { meta, tsv -> [meta, "phix", tsv] }.mix(
        DECONTAMINATE_ASSEMBLIES.out.human_contaminated_contigs_tsv.map { meta, tsv -> [meta, "human", tsv] },
        DECONTAMINATE_ASSEMBLIES.out.host_contaminated_contigs_tsv.map { meta, tsv -> [meta, "contaminant", tsv] }
    )

    /***************************************************************************/
    /*  Cleaned assemblies that fail the following rule:                       */
    /*  - Less than params.short_reads_contig_threshold (default is 2) contigs */
    /***************************************************************************/

    cleaned_contigs.map { meta, assembly_fasta ->
            [meta , ["contigs_count": assembly_fasta.countFasta()], assembly_fasta]
        }
        .branch { _meta, meta2, _assembly_fasta ->
            qc_failed: meta2.contigs_count < params.short_reads_contig_threshold
            qc_passed: meta2.contigs_count >= params.short_reads_contig_threshold
        }
        .set { qc_filtered_assemblies }

    passed_cleaned_contigs = qc_filtered_assemblies.qc_passed.map {
        meta, _meta2, contigs ->
        [ meta, contigs ]
    }

    qc_failed_assemblies = qc_filtered_assemblies.qc_failed.map {
        meta, _meta2, contigs ->
        [meta + ["too_few_contigs": true], contigs]
    }

    emit:
    passed_cleaned_contigs = passed_cleaned_contigs                                         // [ val(meta), path(assembly_fasta) ]
    qc_failed_assemblies   = qc_failed_assemblies                                           // [ val(meta), path(assembly_fasta) ]
    decontaminated_contigs = cleaned_contigs                                                // [ val(meta), path(assembly_fasta) ], all runs
    contaminated_contigs   = contaminated_contigs_by_step                                   // [ val(meta), val(step), path(tsv) ]
    decontamination        = contaminated_contigs_by_step.map { meta, _step, tsv -> [meta, tsv] } // [ val(meta), path(tsv) ]
    decontamination_mqc    = DECONTAMINATE_ASSEMBLIES.out.phix_contaminated_contigs_tsv_mqc.mix(
                                 DECONTAMINATE_ASSEMBLIES.out.human_contaminated_contigs_tsv_mqc,
                                 DECONTAMINATE_ASSEMBLIES.out.host_contaminated_contigs_tsv_mqc
                             )                                                              // [ val(meta), path(mqc_tsv) ]
    versions               = ch_versions
}
