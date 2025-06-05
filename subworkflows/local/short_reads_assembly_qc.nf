include { DECONTAMINATE_CONTIGS as HUMAN_PHIX_DECONTAMINATE_CONTIGS } from '../../subworkflows/ebi-metagenomics/decontaminate_contigs/main'
include { DECONTAMINATE_CONTIGS as HOST_DECONTAMINATE_CONTIGS       } from '../../subworkflows/ebi-metagenomics/decontaminate_contigs/main'
include { SEQKIT_SEQ                                                } from '../../modules/nf-core/seqkit/seq/main'

process PUBLISH_CLEANED_CONTIGS {

    input:
    tuple val(meta), path(cleaned_contigs)

    output:
    tuple val(meta), path("${meta.id}_cleaned.contigs.fa.gz")

    script:
    """
    cp ${cleaned_contigs} ${meta.id}_cleaned.contigs.fa.gz
    """
}

workflow SHORT_READS_ASSEMBLY_QC {

    take:
    assembly               // [ val(meta), path(assembly_fasta) ]

    main:

    ch_versions = Channel.empty()

    /* Len filter using the parameter "short_reads_min_contig_length" */
    SEQKIT_SEQ(
        assembly
    )
    ch_versions = ch_versions.mix(SEQKIT_SEQ.out.versions)

    if ( params.remove_human_phix ) {

        ch_human_phix_refs = Channel
            .fromPath(
                "${params.reference_genomes_folder}/${params.human_phix_ref}",
                checkIfExists: true
            )
            .collect() // TODO Check if collect() is needed here to ensure ref is a val channel
            .map {
                file -> [ ["id": params.human_phix_ref], file ]
            }

        HUMAN_PHIX_DECONTAMINATE_CONTIGS(
            SEQKIT_SEQ.out.fastx,
            ch_human_phix_refs
        )

        prefiltered_assemblies = HUMAN_PHIX_DECONTAMINATE_CONTIGS.out.cleaned_contigs

        ch_versions = ch_versions.mix(HUMAN_PHIX_DECONTAMINATE_CONTIGS.out.versions)

    } else {
        prefiltered_assemblies = SEQKIT_SEQ.out.fastx
    }

    prefiltered_assemblies.branch { meta, assembly_fasta ->
        run_decontamination: meta.contaminant_reference != null
        skip_decontamination: meta.contaminant_reference == null
        }
        .set { subdivided_assemblies }

    subdivided_assemblies.run_decontamination
        .multiMap { meta, assembly_fasta ->
        metagenome: [meta, assembly_fasta]
        contaminant_reference: [meta, meta.contaminant_reference]
        }
        .set { ch_decontamination_input }

    HOST_DECONTAMINATE_CONTIGS(
        ch_decontamination_input.metagenome,
        ch_decontamination_input.contaminant_reference
    )

    cleaned_contigs = subdivided_assemblies.skip_decontamination.mix(HOST_DECONTAMINATE_CONTIGS.out.cleaned_contigs)

    ch_versions = ch_versions.mix(HOST_DECONTAMINATE_CONTIGS.out.versions)

    /***************************************************************************/
    /*  Cleaned assemblies that fail the following rule:                       */
    /*  - Less than params.short_reads_contig_threshold (default is 2) contigs */
    /***************************************************************************/

    cleaned_contigs.map { meta, assembly_fasta -> {
            [meta , ["contigs_count": assembly_fasta.countFasta()], assembly_fasta]
            }
        }
        .branch { _meta, meta2, _assembly_fasta ->
            qc_failed: meta2.contigs_count < params.short_reads_contig_threshold
            qc_passed: meta2.contigs_count >= params.short_reads_contig_threshold
        }
    .set { qc_filtered_assemblies }

    passed_cleaned_contigs = qc_filtered_assemblies.qc_passed.map { meta, _meta2, assembly ->
            [ meta, assembly ]
    }

    qc_failed_assemblies = qc_filtered_assemblies.qc_failed.map { meta, _meta2, assembly ->
        [meta + ["too_few_contigs": true], assembly]
    }

    PUBLISH_CLEANED_CONTIGS(
        passed_cleaned_contigs
    )

    emit:
    passed_cleaned_contigs = passed_cleaned_contigs   // [ val(meta), path(assembly_fasta) ]
    qc_failed_assemblies   = qc_failed_assemblies     // [ val(meta), path(assembly_fasta) ]
    versions               = ch_versions
}
