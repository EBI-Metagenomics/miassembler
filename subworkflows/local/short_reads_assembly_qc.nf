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
    assembly                    // [ val(meta), path(assembly_fasta) ]
    host_reference_genome       // [ val(meta2), path(host_reference_genome) ] | meta2 contains the name of the reference genome

    main:

    ch_versions = Channel.empty()

    /* Len filter using the parameter "short_reads_min_contig_length" */
    SEQKIT_SEQ(
        assembly
    )

    ch_versions = ch_versions.mix(SEQKIT_SEQ.out.versions)

    filtered_contigs = SEQKIT_SEQ.out.fastx

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

    if ( host_reference_genome != null ) {

        ch_blast_host_refs = Channel.fromPath( "${params.blast_reference_genomes_folder}/${host_reference_genome}*", checkIfExists: true)
            .collect().map {
                files -> [ ["id": host_reference_genome], files ]
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

    PUBLISH_CLEANED_CONTIGS(
        filtered_contigs
    )

    emit:
    filtered_contigs = filtered_contigs
    versions         = ch_versions
}
