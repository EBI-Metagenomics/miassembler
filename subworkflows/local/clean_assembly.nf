include { BLAST_BLASTN as BLAST_BLASTN_HUMAN_PHIX } from '../../modules/nf-core/blast/blastn/main'
include { BLAST_BLASTN as BLAST_BLASTN_HOST       } from '../../modules/nf-core/blast/blastn/main'
include { SEQKIT_GREP                             } from '../../modules/nf-core/seqkit/grep/main'
include { SEQKIT_SEQ                              } from '../../modules/nf-core/seqkit/seq/main'

workflow CLEAN_ASSEMBLY {

    take:
    assembly                    // [ val(meta), path(assembly_fasta) ]
    human_phix_reference_genome // [ val(meta2), path(reference_genome) ] | meta2 contains the name of the reference genome
    host_reference_genome       // [ val(meta2), path(reference_genome) ] | meta2 contains the name of the reference genome

    main:

    ch_versions = Channel.empty()

    /* Len filter using the parameter "min_contig_length" */
    SEQKIT_SEQ(
        assembly
    )

    ch_versions = ch_versions.mix(SEQKIT_SEQ.out.versions)

    BLAST_BLASTN_HUMAN_PHIX(
        SEQKIT_SEQ.out.fastx,
        human_phix_reference_genome
    )

    ch_versions = ch_versions.mix(BLAST_BLASTN_HUMAN_PHIX.out.versions.first())

    contaminated_contigs = Channel.empty()
    if ( host_reference_genome != null) {
        ch_blast_host_refs = Channel.fromPath( "$params.blast_reference_genomes_folder/" + host_reference_genome + "*", 
        checkIfExists: true).collect().map {
            files -> [ ["id": host_reference_genome], files ]
        }

        BLAST_BLASTN_HOST(
            SEQKIT_SEQ.out.fastx,
            ch_blast_host_refs
        )

        ch_versions = ch_versions.mix(BLAST_BLASTN_HOST.out.versions.first())

        contaminated_contigs = Channel.of( BLAST_BLASTN_HUMAN_PHIX.out.txt, BLAST_BLASTN_HOST.out.txt )
            .collectFile(storeDir: "${params.outdir}/decontamination/contaminated_contigs.txt", newLine: true)
    } else {
        contaminated_contigs = BLAST_BLASTN_HUMAN_PHIX.out.txt
    }

    SEQKIT_GREP(
        SEQKIT_SEQ.out.fastx,
        contaminated_contigs.map { meta, hits_txt -> {hits_txt }}
    )

    ch_versions = ch_versions.mix(SEQKIT_GREP.out.versions)

    emit:
    filtered_contigs = SEQKIT_GREP.out.filter
    versions         = ch_versions
}
