include { FASTP                                       } from '../../modules/nf-core/fastp/main'
include { READS_BWAMEM2_DECONTAMINATION as PHIX_HUMAN } from '../ebi-metagenomics/reads_bwamem2_decontamination/main'
include { READS_BWAMEM2_DECONTAMINATION as HOST       } from '../ebi-metagenomics/reads_bwamem2_decontamination/main'

workflow PRE_ASSEMBLY_QC {

    take: 
        reads
        humanPhiX_ref_genomes
        host_ref_genomes

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
            humanPhiX_ref_genomes
        )
        
        ch_versions = ch_versions.mix(PHIX_HUMAN.out.versions)
        
        decontaminated_reads = PHIX_HUMAN.out.decontaminated_reads
        if ( host_ref_genomes != null ) {
            ch_bwa_host_refs = Channel.fromPath( "$params.bwa_reference_genomes_folder/" + host_ref_genomes + "*", 
            checkIfExists: true).collect().map {
                files -> [ ["id": host_ref_genomes], files ]
            }
            HOST(
                PHIX_HUMAN.out.decontaminated_reads, 
                ch_bwa_host_refs
            )
            
            ch_versions = ch_versions.mix(HOST.out.versions)
            decontaminated_reads = HOST.out.decontaminated_reads
        }

    emit:
        cleaned_reads = decontaminated_reads
        versions      = ch_versions
}
