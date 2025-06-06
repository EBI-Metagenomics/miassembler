include { DECONTAMINATE_CONTIGS as HUMAN_PHIX_DECONTAMINATE_CONTIGS } from '../../subworkflows/ebi-metagenomics/decontaminate_contigs/main'
include { DECONTAMINATE_CONTIGS as HOST_DECONTAMINATE_CONTIGS       } from '../../subworkflows/ebi-metagenomics/decontaminate_contigs/main'
include { SEQKIT_SEQ                                                } from '../../modules/nf-core/seqkit/seq/main'

workflow SHORT_READS_ASSEMBLY_QC {

    take:
    assembly               // [ val(meta), path(assembly_fasta) ]

    main:

    ch_versions = Channel.empty()

    /* Len filter using the parameter "short_reads_min_contig_length" */
    // TODO add conditional to run filtering only on short reads assemblies
    SEQKIT_SEQ(assembly)
    ch_versions = ch_versions.mix(SEQKIT_SEQ.out.versions)

    if ( params.remove_human_phix ) {

        ch_human_phix_decontamination_input = SEQKIT_SEQ.out.fastx.map {
            // TODO add conditional to choose human+phiX or human decontamination
            meta, assembly_file -> [ [meta, assembly_file ], params.human_phix_ref ]
        }

        HUMAN_PHIX_DECONTAMINATE_CONTIGS(ch_human_phix_decontamination_input)
        prefiltered_assemblies = HUMAN_PHIX_DECONTAMINATE_CONTIGS.out.cleaned_contigs
        ch_versions = ch_versions.mix(HUMAN_PHIX_DECONTAMINATE_CONTIGS.out.versions)

    } else {
        prefiltered_assemblies = SEQKIT_SEQ.out.fastx
    }

    prefiltered_assemblies
        .branch { meta, _contigs ->
            run_decontamination: meta.contaminant_reference != null
            skip_decontamination: meta.contaminant_reference == null
        }
        .set { subdivided_assemblies }

    subdivided_assemblies.run_decontamination
        .map { meta, contigs ->
            [ [meta, contigs], meta.contaminant_reference ]
        }
        .set { ch_decontamination_input }

    HOST_DECONTAMINATE_CONTIGS(ch_decontamination_input)

    cleaned_contigs = subdivided_assemblies.skip_decontamination.mix(
        HOST_DECONTAMINATE_CONTIGS.out.cleaned_contigs
    )

    ch_versions = ch_versions.mix(HOST_DECONTAMINATE_CONTIGS.out.versions)

    //  TODO add check for empty cleaned contigs to SR and/or LR workflow

    emit:
    cleaned_contigs = cleaned_contigs  // [ val(meta), path(assembly_fasta) ]
    versions        = ch_versions
}
