include { BLAST_BLASTN } from '../modules/nf-core/blast/blastn/main'
include { SEQKIT_GREP } from '../modules/nf-core/seqkit/grep/main'

workflow CLEAN_ASSEMBLY {

    take:
    assembly         // [ val(meta), path(assembly_fasta) ]
    reference_genome // [ val(meta), path(reference_genome) ]

    main:

    ch_versions = Channel.empty()

    BLAST_BLASTN(
        assembly,
        reference_genome
    )

    ch_versions = ch_versions.mix(BLAST_BLASTN.out.versions.first())

    SEQKIT_GREP(
        assembly.join( BLAST_BLASTN.out.txt )
    )

    ch_versions = ch_versions.mix(SEQKIT_GREP.out.versions.first())

    output:
    filtered_contigs = SEQKIT_GREP.out.filter
    versions         = ch_versions
}
