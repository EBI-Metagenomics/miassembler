include { PROKKA         } from '../../modules/local/prokka'
include { DIAMOND_BLASTP } from '../../modules/nf-core/diamond/blastp/main'
include { IDEEL          } from '../../modules/local/ideel'
include { PROOVFRAME_FIX } from '../../modules/nf-core/proovframe/fix/main'
include { PROOVFRAME_MAP } from '../../modules/nf-core/proovframe/map/main'

workflow FRAMESHIFT_CORRECTION {
    take:
    contigs                   // [ val(meta), path(contigs) ]

    main:
    ch_versions = Channel.empty()
    contigs.view()

    PROKKA(
      contigs
    )
    ch_versions = ch_versions.mix(PROKKA.out.versions)

    DIAMOND_BLASTP(
      PROKKA.out.prokka_faa,
      "/home/germanab/Downloads/uniprot_sprot.fasta.gz",
      "txt",
      ""
    )
    ch_versions = ch_versions.mix(DIAMOND_BLASTP.out.versions)

    // PROOVFRAME_MAP(
    //     reads
    // )
    // ch_versions = ch_versions.mix(PROOVFRAME_MAP.out.versions)

    emit:
    // corrected_contigs  = DIAMOND.out.tsv
    corrected_contigs = PROKKA.out.prokka_faa
    versions = ch_versions
}
