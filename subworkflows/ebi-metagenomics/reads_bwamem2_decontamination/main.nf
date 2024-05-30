include { BWAMEM2_MEM as BWAMEM2_MEM_DECONT } from '../../../modules/ebi-metagenomics/bwamem2/mem/main'
include { SAMTOOLS_BAM2FQ                   } from '../../../modules/ebi-metagenomics/samtools/bam2fq/main'


workflow READS_BWAMEM2_DECONTAMINATION {

    take:
    ch_reads       // channel (mandatory): [ val(meta), [ path(reads) ] ]
    ch_reference   // channel (mandatory): [ val(meta2), path(ref_index) ] | meta2 contains the name of the reference genome


    main:

    ch_versions = Channel.empty()

    BWAMEM2_MEM(ch_reads, ch_reference)
    ch_versions = ch_versions.mix(BWAMEM2_MEM.out.versions.first())

    SAMTOOLS_BAM2FQ( BWAMEM2_MEM.out.bam.map { meta, bam, _ -> [ meta, bam, meta.single_end == false ] } )
    ch_versions = ch_versions.mix(SAMTOOLS_BAM2FQ.out.versions.first())

    emit:
    decontaminated_reads = SAMTOOLS_BAM2FQ.out.reads  // channel: [ val(meta), [ path(decont_reads) ]]
    versions = ch_versions                            // channel: [ versions.yml ]

}
