include { FASTP                                              } from '../../modules/nf-core/fastp/main'
include { BWAMEM2DECONTNOBAMS as HUMAN_READS_DECONTAMINATION } from '../../modules/ebi-metagenomics/bwamem2decontnobams/main'
include { BWAMEM2DECONTNOBAMS as PHIX_READS_DECONTAMINATION  } from '../../modules/ebi-metagenomics/bwamem2decontnobams/main'
include { BWAMEM2DECONTNOBAMS as HOST_READS_DECONTAMINATION  } from '../../modules/ebi-metagenomics/bwamem2decontnobams/main'

workflow SHORT_READS_QC {

    take:
    reads            // [ val(meta), path(reads) ]

    main:
    ch_versions = Channel.empty()

    FASTP(
        reads,
        [],
        false,
        false,
        false,
        false
    )

    ch_versions = ch_versions.mix(FASTP.out.versions)

    /***************************************************************************/
    /* Perform decontamination from human sequences if requested               */
    /***************************************************************************/
    FASTP.out.reads
        .branch { meta, _reads ->
            run_decontamination: meta.human_reference != null
            skip_decontamination: meta.human_reference == null
        }
        .set { human_subdivided_reads }

    human_subdivided_reads.run_decontamination
        .multiMap { meta, reads ->
            reads: [meta, reads]
            reference: [ [id:"human"], file("${params.reference_genomes_folder}/${meta.human_reference}.*") ]
        }
        .set { ch_human_decontamination_input }

    HUMAN_READS_DECONTAMINATION(
            ch_human_decontamination_input.reads,
            ch_human_decontamination_input.reference
        )

    ch_versions = ch_versions.mix(HUMAN_READS_DECONTAMINATION.out.versions)

    human_cleaned_reads = human_subdivided_reads.skip_decontamination.mix(
        HUMAN_READS_DECONTAMINATION.out.decont_reads
    )

    /***************************************************************************/
    /* Perform decontamination from PhiX sequences if requested               */
    /***************************************************************************/
    human_cleaned_reads
        .branch { meta, _reads ->
            run_decontamination: meta.phix_reference != null
            skip_decontamination: meta.phix_reference == null
        }
        .set { phix_subdivided_reads }

    phix_subdivided_reads.run_decontamination
        .multiMap { meta, reads ->
            reads: [meta, reads]
            reference: [ [id:"phix"], file("${params.reference_genomes_folder}/${meta.phix_reference}.*") ]
        }
        .set { ch_phix_decontamination_input }

    PHIX_READS_DECONTAMINATION(
            ch_phix_decontamination_input.reads,
            ch_phix_decontamination_input.reference
        )

    ch_versions = ch_versions.mix(PHIX_READS_DECONTAMINATION.out.versions)

    phix_cleaned_reads = phix_subdivided_reads.skip_decontamination.mix(
        PHIX_READS_DECONTAMINATION.out.decont_reads
    )

    /***************************************************************************/
    /* Perform decontamination from arbitrary contaminant sequences            */
    /***************************************************************************/
    phix_cleaned_reads
        .branch { meta, _reads ->
            run_decontamination: meta.contaminant_reference != null
            skip_decontamination: meta.contaminant_reference == null
        }
        .set { subdivided_reads }

    subdivided_reads.run_decontamination
        .multiMap { meta, reads ->
            reads: [meta, reads]
            reference: [ [id:file(meta.contaminant_reference).baseName], file("${params.reference_genomes_folder}/${meta.contaminant_reference}.*") ]
        }
        .set { ch_decontamination_input }

    HOST_READS_DECONTAMINATION(
            ch_decontamination_input.reads,
            ch_decontamination_input.reference
        )

    ch_versions = ch_versions.mix(HOST_READS_DECONTAMINATION.out.versions)

    decontaminated_reads = subdivided_reads.skip_decontamination.mix(
        HOST_READS_DECONTAMINATION.out.decont_reads
    )

    emit:
    qc_reads   = decontaminated_reads
    fastp_json = FASTP.out.json
    versions   = ch_versions
}
