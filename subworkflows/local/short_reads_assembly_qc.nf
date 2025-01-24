include { BLAST_BLASTN as BLAST_BLASTN_HUMAN_PHIX } from '../../modules/nf-core/blast/blastn/main'
include { BLAST_BLASTN as BLAST_BLASTN_HOST       } from '../../modules/nf-core/blast/blastn/main'
include { SEQKIT_GREP as SEQKIT_GREP_HUMAN_PHIX   } from '../../modules/nf-core/seqkit/grep/main'
include { SEQKIT_GREP as SEQKIT_GREP_HOST         } from '../../modules/nf-core/seqkit/grep/main'
include { SEQKIT_SEQ                              } from '../../modules/nf-core/seqkit/seq/main'

process PUBLISH_CLEANED_CONTIGS {

    input:
    tuple val(meta), path(cleaned_contigs)

    output:
    tuple val(meta), path("${meta.id}_cleaned.contigs.fa.gz")

    script:
    """
    cp ${cleaned_contigs} ${meta.id}_cleaned.contigs.fa.gz
    """
}

workflow SHORT_READS_ASSEMBLY_QC {

    take:
    assembly               // [ val(meta), path(assembly_fasta) ]
    reference_genome       // [ val(meta2), path(reference_genome) ] | meta2 contains the name of the reference genome

    main:

    def ch_versions = Channel.empty()
    def ch_blast_human_phix_refs = Channel.empty()
    def ch_blast_host_refs = Channel.empty()

    /* Len filter using the parameter "short_reads_min_contig_length" */
    SEQKIT_SEQ(
        assembly
    )

    ch_versions = ch_versions.mix(SEQKIT_SEQ.out.versions)

    def filtered_contigs = SEQKIT_SEQ.out.fastx

    if ( params.remove_human_phix ) {

        ch_blast_human_phix_refs = Channel.fromPath( "${params.blast_reference_genomes_folder}/${params.human_phix_blast_index_name}*", checkIfExists: true)
            .collect().map {
                files -> [ ["id": params.human_phix_blast_index_name], files ]
            }

        BLAST_BLASTN_HUMAN_PHIX(
            SEQKIT_SEQ.out.fastx,
            ch_blast_human_phix_refs
        )

        ch_versions = ch_versions.mix(BLAST_BLASTN_HUMAN_PHIX.out.versions.first())

        SEQKIT_GREP_HUMAN_PHIX(
            filtered_contigs.join( BLAST_BLASTN_HUMAN_PHIX.out.txt )
        )

        filtered_contigs = SEQKIT_GREP_HUMAN_PHIX.out.filter

        ch_versions = ch_versions.mix(SEQKIT_GREP_HUMAN_PHIX.out.versions)
    }
    
    // The cleaned contigs are those that have been filtered, but they will be further cleaned if a reference genome is set.
    def cleaned_contigs = filtered_contigs

    if ( reference_genome != null ) {

        ch_blast_host_refs = Channel.fromPath( "${params.blast_reference_genomes_folder}/${reference_genome}*", checkIfExists: true)
            .collect().map {
                files -> [ ["id": reference_genome], files ]
            }

        BLAST_BLASTN_HOST(
            filtered_contigs,
            ch_blast_host_refs
        )

        ch_versions = ch_versions.mix(BLAST_BLASTN_HOST.out.versions.first())

        SEQKIT_GREP_HOST(
            filtered_contigs.join( BLAST_BLASTN_HOST.out.txt )
        )

        cleaned_contigs = SEQKIT_GREP_HOST.out.filter

        ch_versions = ch_versions.mix(SEQKIT_GREP_HOST.out.versions)
    }

    /***************************************************************************/
    /*  Cleaned assemblies that fail the following rule:                       */
    /*  - Less than params.short_reads_contig_threshold (default is 2) contigs */
    /***************************************************************************/

    cleaned_contigs.map { meta, assembly_fasta -> {
           [meta , ["contigs_count": assembly_fasta.countFasta()], assembly_fasta]
            }
        }
        .branch { meta, meta2, assembly_fasta ->
            qc_failed: meta2.contigs_count < params.short_reads_contig_threshold
            qc_passed: meta2.contigs_count >= params.short_reads_contig_threshold
        }
    .set { qc_filtered_assemblies }
    
    passed_cleaned_contigs = qc_filtered_assemblies.qc_passed.map { meta, _meta2, assembly -> 
            [ meta, assembly ]
    }

    qc_failed_assemblies = qc_filtered_assemblies.qc_failed.map { meta, _meta2, assembly -> 
        [meta + ["too_few_contigs": true], assembly] 
    }

    PUBLISH_CLEANED_CONTIGS(
        passed_cleaned_contigs
    )

    emit:
    passed_cleaned_contigs = passed_cleaned_contigs // tuple(meta)
    qc_failed_assemblies   = qc_failed_assemblies // tuple(meta)
    versions               = ch_versions
}
