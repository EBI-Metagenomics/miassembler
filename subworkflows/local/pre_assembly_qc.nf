include { FASTP                                                       } from '../../modules/nf-core/fastp/main'
include { READS_BWAMEM2_DECONTAMINATION as HUMAN_PHIX_DECONTAMINATION } from '../ebi-metagenomics/reads_bwamem2_decontamination/main'
include { READS_BWAMEM2_DECONTAMINATION as HOST_DECONTAMINATION       } from '../ebi-metagenomics/reads_bwamem2_decontamination/main'

workflow PRE_ASSEMBLY_QC {

    take: 
        reads                      // [ val(meta), path(reads) ]
        human_phix_reference_genome
        host_reference_genome      // [ val(meta2), path(reference_genome) ] | meta2 contains the name of the reference genome

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

        HUMAN_PHIX_DECONTAMINATION(
            FASTP.out.reads, 
            human_phix_reference_genome
        )
        
        ch_versions = ch_versions.mix(HUMAN_PHIX_DECONTAMINATION.out.versions)
        
        decontaminated_reads = HUMAN_PHIX_DECONTAMINATION.out.decontaminated_reads
        if ( host_reference_genome != null ) {
            ch_bwa_host_refs = Channel.fromPath( "$params.bwa_reference_genomes_folder/" + host_reference_genome + "*", 
            checkIfExists: true).collect().map {
                files -> [ ["id": host_reference_genome], files ]
            }
            HOST_PHIX_DECONTAMINATION(
                HUMAN_PHIX_DECONTAMINATION.out.decontaminated_reads, 
                ch_bwa_host_refs
            )
            
            ch_versions = ch_versions.mix(HOST_PHIX_DECONTAMINATION.out.versions)
            decontaminated_reads = HOST_PHIX_DECONTAMINATION.out.decontaminated_reads
        }

    emit:
        cleaned_reads = decontaminated_reads
        versions      = ch_versions
}
