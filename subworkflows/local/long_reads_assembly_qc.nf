include { DECONTAMINATE_ASSEMBLIES } from '../../subworkflows/ebi-metagenomics/decontaminate_assemblies/main'

workflow LONG_READS_ASSEMBLY_QC {

    take:
    assembly            // [ val(meta), path(assembly_fasta) ]

    main:
    ch_versions = channel.empty()

    /***************************************************************************/
    /* Remove human and host contigs (each step runs only if requested)        */
    /* PhiX decontamination is not used for long reads                         */
    /***************************************************************************/
    DECONTAMINATE_ASSEMBLIES(
        assembly.map { meta, contigs -> [meta + [phix_reference: null], contigs] }
    )
    ch_versions = ch_versions.mix(DECONTAMINATE_ASSEMBLIES.out.versions)

    def contaminated_contigs_by_step = DECONTAMINATE_ASSEMBLIES.out.human_contaminated_contigs_tsv.map { meta, tsv -> [meta, "human", tsv] }.mix(
        DECONTAMINATE_ASSEMBLIES.out.host_contaminated_contigs_tsv.map { meta, tsv -> [meta, "contaminant", tsv] }
    )

    emit:
    contigs              = DECONTAMINATE_ASSEMBLIES.out.cleaned_contigs                         // [ val(meta), path(assembly_fasta) ]
    contaminated_contigs = contaminated_contigs_by_step                                         // [ val(meta), val(step), path(tsv) ]
    decontamination      = contaminated_contigs_by_step.map { meta, _step, tsv -> [meta, tsv] } // [ val(meta), path(tsv) ]
    decontamination_mqc  = DECONTAMINATE_ASSEMBLIES.out.human_contaminated_contigs_tsv_mqc.mix(
                               DECONTAMINATE_ASSEMBLIES.out.host_contaminated_contigs_tsv_mqc
                           )                                                                    // [ val(meta), path(mqc_tsv) ]
    versions             = ch_versions
}
