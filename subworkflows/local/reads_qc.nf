include { FASTP                                                       } from '../../modules/nf-core/fastp/main'
include { READS_BWAMEM2_DECONTAMINATION as HUMAN_PHIX_DECONTAMINATION } from '../ebi-metagenomics/reads_bwamem2_decontamination/main'
include { READS_BWAMEM2_DECONTAMINATION as HOST_DECONTAMINATION       } from '../ebi-metagenomics/reads_bwamem2_decontamination/main'

workflow READS_QC {

    take:
    reads                 // [ val(meta), path(reads) ]
    host_reference_genome // [ val(meta2), path(reference_genome) ] | meta2 contains the name of the reference genome
    isMetatranscriptomic  // true/false

    main:
    ch_versions = Channel.empty()

    format_reads = { meta, reads ->
        if (reads.size() == 1) {
            return tuple(meta + [single_end: true], reads)
        }
        else {
            return tuple(meta + [single_end: false], reads)
        }
    }

    // TO BE CHANGED ONCE SAMPLESHEETS ARE SET UP
    // formatted_reads = samplesheet.map(format_reads)
    formatted_reads = reads.map(format_reads)

    // paired/single_end should be removed in favour of samplesheets input

    FASTP(
        formatted_reads,
        [],
        false,
        false,
        isMetatranscriptomic
    )

    ch_versions = ch_versions.mix(FASTP.out.versions)

    decontaminated_reads = channel.empty()

    if ( params.remove_human_phix ) {

        ch_bwamem2_human_phix_refs = Channel.fromPath( "$params.bwamem2_reference_genomes_folder/${params.human_phix_blast_index_name}*", checkIfExists: true)
            .collect().map {
                files -> [ ["id": params.human_phix_blast_index_name], files ]
            }

        HUMAN_PHIX_DECONTAMINATION(
            FASTP.out.reads,
            ch_bwamem2_human_phix_refs
        )

        ch_versions = ch_versions.mix(HUMAN_PHIX_DECONTAMINATION.out.versions)

        decontaminated_reads = HUMAN_PHIX_DECONTAMINATION.out.decontaminated_reads

    } else {
        decontaminated_reads = FASTP.out.reads
    }

    if ( host_reference_genome != null ) {

        ch_bwamem2_host_refs = Channel.fromPath( "${params.bwamem2_reference_genomes_folder}/${host_reference_genome}*", checkIfExists: true)
            .collect().map {
                files -> [ ["id": host_reference_genome], files ]
            }

        HOST_DECONTAMINATION(
            HUMAN_PHIX_DECONTAMINATION.out.decontaminated_reads,
            ch_bwamem2_host_refs
        )

        ch_versions = ch_versions.mix(HOST_DECONTAMINATION.out.versions)

        decontaminated_reads = HOST_DECONTAMINATION.out.decontaminated_reads
    }

    emit:
    qc_reads = decontaminated_reads
    versions = ch_versions
}
