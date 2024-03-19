include { FASTP                                       } from '../../modules/nf-core/fastp/main'
include { READS_BWAMEM2_DECONTAMINATION as PHIX_HUMAN } from '../ebi-metagenomics/reads_bwamem2_decontamination/main'
include { READS_BWAMEM2_DECONTAMINATION as HOST       } from '../ebi-metagenomics/reads_bwamem2_decontamination/main'

workflow PRE_ASSEMBLY_QC {

    take: 
        reads
        ref_genomes

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
            false 
        )
        
        ch_versions = ch_versions.mix(FASTP.out.versions)

        PHIX_HUMAN(
            FASTP.out.reads, 
            ref_genomes
        )
        
        ch_versions = ch_versions.mix(PHIX_HUMAN.out.versions)
        decontaminated_reads = Channel.empty()

        if (ref_genomes.count() == 2) {
            HOST(
                PHIX_HUMAN.out.decontaminated_reads, 
                ref_genomes[1]
            )
            
            ch_versions = ch_versions.mix(HOST.out.versions)
            decontaminated_reads = HOST.out.decontaminated_reads
            HOST.out.decontaminated_reads.view()
        }
        else {
            decontaminated_reads = PHIX_HUMAN.out.decontaminated_reads
        }

    emit:
        cleaned_reads = decontaminated_reads
        versions      = ch_versions
}
