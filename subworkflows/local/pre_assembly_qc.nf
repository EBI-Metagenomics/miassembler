include { FASTP } from '../../modules/nf-core/fastp/main
include { READS_QC } from 'reads_bwamem_decont/main'

workflow PRE_ASSEMBLY_QC {

    take:          // reads? // from Kate: [[meta], [reads]]

    main:
    ch_versions = Channel.empty()

    FASTP(         // look at output of fetchtool or samplesheet

    )

    ch_versions = ch_versions.mix(FASTP.out.versions)

    READS_QC(
        //reads               // tuple(meta, reads) // the ones coming from FASTP
        //ref_genome          // path(reference_genome)
        //ref_genome_index    // path(reference_genome_index
    )

    ch_versions = ch_versions.mix(READS_QC.out.versions)

    emit:
    cleaned_reads = READS_QC.out.filter
    versions      = ch_versions
}
