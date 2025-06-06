include { DECONTAMINATE_CONTIGS as HUMAN_DECONTAMINATE_CONTIGS } from '../../subworkflows/ebi-metagenomics/decontaminate_contigs/main'
include { DECONTAMINATE_CONTIGS as HOST_DECONTAMINATE_CONTIGS  } from '../../subworkflows/ebi-metagenomics/decontaminate_contigs/main'

workflow LONG_READS_ASSEMBLY_QC {

    take:
    assembly            // [ val(meta), path(assembly_fasta) ]

    main:
    ch_versions = Channel.empty()

    if ( params.remove_human ) {

        ch_human_decontamination_input = assembly.map {
            // TODO decide if params.human_phix_ref is filename or path
            meta, assembly_file -> [ [meta, assembly_file ], params.human_ref ]
        }

        HUMAN_DECONTAMINATE_CONTIGS(ch_human_decontamination_input)
        prefiltered_assemblies = HUMAN_DECONTAMINATE_CONTIGS.out.cleaned_contigs
        ch_versions = ch_versions.mix(HUMAN_DECONTAMINATE_CONTIGS.out.versions)

    } else {
        prefiltered_assemblies = assembly
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

    emit:
    contigs  = cleaned_contigs
    versions = ch_versions
}