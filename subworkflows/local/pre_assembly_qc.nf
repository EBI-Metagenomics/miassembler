include { FASTP as FASTP_2                            } from '../../modules/nf-core/fastp/main'
include { READS_BWAMEM2_DECONTAMINATION as PHIX_HUMAN } from '../ebi-metagenomics/reads_bwamem2_decontamination/main'
include { READS_BWAMEM2_DECONTAMINATION as HOST       } from '../ebi-metagenomics/reads_bwamem2_decontamination/main'

workflow PRE_ASSEMBLY_QC {

    take: 
        reads
        ref_genomes

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

        // read_input = samplesheet.map(format_reads)
        read_input = format_reads
        ref_genomes.view()

        // FASTP_2( 
        //     read_input, 
        //     [], 
        //     false, 
        //     false 
        // )
        // FASTP.out.reads.view()

        // ch_versions = ch_versions.mix(FASTP.out.versions)

        // PHIX_HUMAN(
        //     FASTP.out.reads, 
        //     ref_genomes[0]
        // )

        // ch_versions = ch_versions.mix(PHIX_HUMAN.out.versions)
        // decontaminated_reads = Channel.empty()

        // if (myChannel.count() == 2) {
        //     HOST(
        //         PHIX_HUMAN.out.reads, 
        //         ref_genomes[1]
        //     )
        //     ch_versions = ch_versions.mix(HOST.out.versions)
        //     decontaminated_reads = HOST.out.reads
        // }
        // else {
        //     decontaminated_reads = PHIX_HUMAN.out.reads
        // }

    emit:
    // cleaned_reads = decontaminated_reads
        cleaned_reads = reads
        versions      = ch_versions
}
