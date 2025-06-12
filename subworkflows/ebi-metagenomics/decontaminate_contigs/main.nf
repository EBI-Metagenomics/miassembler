include { MINIMAP2_ALIGN } from '../../../modules/nf-core/minimap2/align/main'
include { FILTERPAF      } from '../../../modules/ebi-metagenomics/filterpaf/main'
include { SEQKIT_GREP    } from '../../../modules/nf-core/seqkit/grep/main'

workflow DECONTAMINATE_CONTIGS {

    // ================================================================================
    // Decontamination
    //  Remove contigs with query coverage ≥ min_qcov
    //    AND
    //  percentage identity ≥ min_pid
    //
    // TODO: Long PacBio HiFi contigs (8+ Mbp) may have large partial
    // alignments (3+ Mbp) to contaminants, making absolute length more
    // informative than coverage percentage.
    // To address this, we could add additional filtering criteria:
    //  alignment > min_align_len bp
    // ================================================================================

    take:
    contigs_and_reference      // [ [ meta, path(assembly_fasta)], path(reference_fasta) ]

    main:
    ch_versions = Channel.empty()

    contigs_and_reference
        .multiMap { assembly, contaminant_reference ->
            contigs: assembly
            reference: [[id: file(contaminant_reference).baseName], contaminant_reference ]
        }
        .set { minimap2_input_ch }

    MINIMAP2_ALIGN(
        minimap2_input_ch.contigs,
        minimap2_input_ch.reference,
        false,             // bam_format
        false,             // bam_index_extension
        false,             // cigar_paf_format
        false              // cigar_bam
    )
    ch_versions = ch_versions.mix(MINIMAP2_ALIGN.out.versions.first())

    FILTERPAF(
        MINIMAP2_ALIGN.out.paf
    )
    ch_versions = ch_versions.mix(FILTERPAF.out.versions.first())

    minimap2_input_ch.contigs
        .join(FILTERPAF.out.mapped_contigs_txt)
        .multiMap { meta, assembly, mapped_contigs_txt ->
            contigs: [meta, assembly]
            pattern: mapped_contigs_txt
        }
        .set { seqkit_grep_input_ch }

    SEQKIT_GREP(
        seqkit_grep_input_ch.contigs,
        seqkit_grep_input_ch.pattern,
    )
    ch_versions = ch_versions.mix(SEQKIT_GREP.out.versions.first())

    emit:
    cleaned_contigs = SEQKIT_GREP.out.filter
    versions        = ch_versions
}
