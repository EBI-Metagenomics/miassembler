include { PROKKA as PROKKA_BEFORE } from '../../modules/local/prokka'
include { PROKKA as PROKKA_AFTER  } from '../../modules/local/prokka'
include { DIAMOND_MAKEDB        } from '../../modules/nf-core/diamond/makedb/main'
include { DIAMOND_BLASTP as DIAMOND_BLASTP_BEFORE } from '../../modules/nf-core/diamond/blastp/main'
include { DIAMOND_BLASTP as DIAMOND_BLASTP_AFTER } from '../../modules/nf-core/diamond/blastp/main'
include { IDEEL as IDEEL_BEFORE } from '../../modules/local/ideel'
include { IDEEL as IDEEL_AFTER  } from '../../modules/local/ideel'
include { PROOVFRAME_FIX        } from '../../modules/nf-core/proovframe/fix/main'
include { PROOVFRAME_MAP        } from '../../modules/nf-core/proovframe/map/main'

workflow FRAMESHIFT_CORRECTION {
    take:
    contigs                   // [ val(meta), path(contigs) ]

    main:
    ch_versions = Channel.empty()

    PROKKA_BEFORE(
        contigs
    )
    ch_versions = ch_versions.mix(PROKKA_BEFORE.out.versions)

    DIAMOND_MAKEDB(
        [[ id:'diamond_db' ], file(params.diamond_db, checkIfExists: true) ], [], [], []
    )
    ch_versions = ch_versions.mix(DIAMOND_MAKEDB.out.versions)

    DIAMOND_BLASTP_BEFORE(
      PROKKA_BEFORE.out.prokka_faa,
      DIAMOND_MAKEDB.out.db,
      "txt",
      "qlen slen"
    )
    ch_versions = ch_versions.mix(DIAMOND_BLASTP_BEFORE.out.versions)

    IDEEL_BEFORE(
        DIAMOND_BLASTP_BEFORE.out.txt
    )

    PROOVFRAME_MAP(
        PROKKA_BEFORE.out.prokka_faa,
        contigs
    )
    ch_versions = ch_versions.mix(PROOVFRAME_MAP.out.versions)

    PROOVFRAME_FIX(
        contigs, 
        PROOVFRAME_MAP.out.tsv
    )
    ch_versions = ch_versions.mix(PROOVFRAME_FIX.out.versions)

    PROKKA_AFTER(
        PROOVFRAME_FIX.out.out_fa
    )
    ch_versions = ch_versions.mix(PROKKA_AFTER.out.versions)

    DIAMOND_BLASTP_AFTER(
      PROKKA_AFTER.out.prokka_faa,
      DIAMOND_MAKEDB.out.db,
      "txt",
      "qlen slen"
    )
    ch_versions = ch_versions.mix(DIAMOND_BLASTP_AFTER.out.versions)

    IDEEL_AFTER(
        DIAMOND_BLASTP_AFTER.out.txt
    )

    emit:
    corrected_contigs  = PROOVFRAME_FIX.out.out_fa
    versions = ch_versions
}
