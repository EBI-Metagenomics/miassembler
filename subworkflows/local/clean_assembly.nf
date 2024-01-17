include { BLAST_BLASTN } from '../../modules/nf-core/blast/blastn/main'
include { SEQKIT_GREP } from '../../modules/nf-core/seqkit/grep/main'
include { SEQKIT_SEQ } from '../../modules/nf-core/seqkit/seq/main'

workflow CLEAN_ASSEMBLY {

    take:
    assembly         // [ val(meta), path(assembly_fasta) ]
    reference_genome // [ val(meta2), path(reference_genome) ] | meta2 contains the name of the reference genome

    main:

    ch_versions = Channel.empty()

    /* Len filter using the parameter "min_contig_length" */
    SEQKIT_SEQ(
        assembly
    )

    BLAST_BLASTN(
        SEQKIT_SEQ.out.fastx,
        reference_genome
    )

    ch_versions = ch_versions.mix(BLAST_BLASTN.out.versions.first())

    SEQKIT_GREP(
        SEQKIT_SEQ.out.fastx,
        BLAST_BLASTN.out.txt.map { meta, hits_txt -> hits_txt }
    )

    ch_versions = ch_versions.mix(SEQKIT_GREP.out.versions.first())

    emit:
    filtered_contigs = SEQKIT_GREP.out.filter
    versions         = ch_versions
}
