include { BLAST_BLASTN as BLAST_BLASTN_HUMAN_PHIX } from '../../modules/nf-core/blast/blastn/main'
include { BLAST_BLASTN as BLAST_BLASTN_HOST       } from '../../modules/nf-core/blast/blastn/main'
include { SEQKIT_GREP                             } from '../../modules/nf-core/seqkit/grep/main'
include { SEQKIT_SEQ                              } from '../../modules/nf-core/seqkit/seq/main'
include { PUBLISH_FILE                            } from '../../modules/local/publish_file'

workflow ASSEMBLY_QC {

    take:
    assembly                    // [ val(meta), path(assembly_fasta) ]
    host_reference_genome       // [ val(meta2), path(reference_genome) ] | meta2 contains the name of the reference genome

    main:

    ch_versions = Channel.empty()

    /* Len filter using the parameter "min_contig_length" */
    SEQKIT_SEQ(
        assembly
    )

    ch_versions = ch_versions.mix(SEQKIT_SEQ.out.versions)

    contaminated_contigs = channel.empty()

    if ( params.remove_human_phix ) {

        ch_blast_human_phix_refs = Channel.fromPath( "$params.blast_reference_genomes_folder/${params.human_phix_blast_index_name}*", checkIfExists: true)
            .collect().map {
                files -> [ ["id": params.human_phix_blast_index_name], files ]
            }

        BLAST_BLASTN_HUMAN_PHIX(
            SEQKIT_SEQ.out.fastx,
            ch_blast_human_phix_refs
        )

        ch_versions = ch_versions.mix(BLAST_BLASTN_HUMAN_PHIX.out.versions.first())

        contaminated_contigs = BLAST_BLASTN_HUMAN_PHIX.out.txt
    }

    if ( host_reference_genome != null) {

        ch_blast_host_refs = Channel.fromPath( "${params.blast_reference_genomes_folder}/${host_reference_genome}*", checkIfExists: true)
            .collect().map {
                files -> [ ["id": host_reference_genome], files ]
            }

        BLAST_BLASTN_HOST(
            SEQKIT_SEQ.out.fastx,
            ch_blast_host_refs
        )

        ch_versions = ch_versions.mix(BLAST_BLASTN_HOST.out.versions.first())

        contaminated_contigs = Channel.of( BLAST_BLASTN_HUMAN_PHIX.out.txt, BLAST_BLASTN_HOST.out.txt )
            .collectFile(name: "contaminated_contigs_host.txt", newLine: true)
    } else {
        contaminated_contigs = BLAST_BLASTN_HUMAN_PHIX.out.txt
    }

    collected_contigs = contaminated_contigs.map { meta, hits_txt -> hits_txt }.collectFile(name: "decontaminated.txt", newLine: true)
    output_path = contaminated_contigs.map { meta, _ -> "assembly/$meta.assembler/$meta.assembler_version/decontamination/"}
    PUBLISH_FILE(collected_contigs, output_path)

    SEQKIT_GREP(
        SEQKIT_SEQ.out.fastx,
        contaminated_contigs.map { meta, hits_txt -> { hits_txt }}
    )

    ch_versions = ch_versions.mix(SEQKIT_GREP.out.versions)

    emit:
    filtered_contigs = SEQKIT_GREP.out.filter
    versions         = ch_versions
}
