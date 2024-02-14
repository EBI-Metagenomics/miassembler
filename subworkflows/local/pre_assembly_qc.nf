include { FASTP               } from '../../modules/nf-core/fastp/main'
include { READS_BWAMEM_DECONT } from 'reads_bwamem_decont/main'

workflow PRE_ASSEMBLY_QC {

    take: 
        reads
        ref_genome
        ref_genome_index

    main:
        ch_versions = Channel.empty()

    format_reads = { meta, fq1, fq2 ->
        if (fq2 == []) {
            return tuple(meta, [fq1])
        }
        else {
            return tuple(meta, [fq1, fq2])
        }
    }

    read_input = samplesheet.map(format_reads)

    FASTP(
        read_input,
        [],
        false,
        false
    )

    ch_versions = ch_versions.mix(FASTP.out.versions)

    DECONTAMINATION(
        reads,
        ref_genome,
        ref_genome_index,
    )

    ch_versions = ch_versions.mix(DECONTAMINATION.out.versions)

    emit:
    cleaned_reads = DECONTAMINATION.out.filter
    versions      = ch_versions
}
