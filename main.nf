#!/usr/bin/env nextflow
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ebi-metagenomics/miassembler
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Github : https://github.com/ebi-metagenomics/miassembler
----------------------------------------------------------------------------------------
*/

include { validateParameters; paramsHelp } from 'plugin/nf-schema'

include { MIASSEMBLER } from './workflows/miassembler'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    OUTPUT PATHS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// <outdir>/PRJEB12/PRJEB12345
def studyDir(Map meta) {
    "${meta.study_accession.substring(0, 7)}/${meta.study_accession}"
}

// <outdir>/PRJEB12/PRJEB12345/SRR1234/SRR1234567
def runDir(Map meta) {
    "${studyDir(meta)}/${meta.id.substring(0, 7)}/${meta.id}"
}

// <outdir>/PRJEB12/PRJEB12345/SRR1234/SRR1234567/assembly/<assembler>/<assembler_version>
def assemblyDir(Map meta) {
    "${runDir(meta)}/assembly/${meta.assembler}/${meta.assembler_version}"
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN THE PIPELINE
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow {

    main:

    if (params.help) {
        log.info paramsHelp(
            command: "nextflow run ebi-metagenomics/miassembler --help",
            showHidden: params.show_hidden,
        )
        exit 0
    }

    validateParameters()

    MIASSEMBLER()

    publish:
    fastqc_before    = MIASSEMBLER.out.fastqc_before
    fastqc_after     = MIASSEMBLER.out.fastqc_after
    fastp            = MIASSEMBLER.out.fastp
    fastplong        = MIASSEMBLER.out.fastplong
    assembler_params = MIASSEMBLER.out.assembler_params
    contigs          = MIASSEMBLER.out.contigs
    decontamination  = MIASSEMBLER.out.decontamination
    coverage         = MIASSEMBLER.out.coverage
    quast            = MIASSEMBLER.out.quast
    multiqc_study    = MIASSEMBLER.out.multiqc_study
    multiqc_run      = MIASSEMBLER.out.multiqc_run
    versions         = MIASSEMBLER.out.versions
    assembled_runs   = MIASSEMBLER.out.assembled_runs
    qc_failed_runs   = MIASSEMBLER.out.qc_failed_runs
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    PUBLISHED OUTPUTS (relative to params.outdir)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

output {

    // Reads QC //
    fastqc_before {
        path { meta, _files -> "${runDir(meta)}/qc/fastqc/before" }
    }
    fastqc_after {
        path { meta, _files -> "${runDir(meta)}/qc/fastqc/after" }
    }
    fastp {
        path { meta, _files -> "${runDir(meta)}/qc/fastp" }
    }
    fastplong {
        path { meta, _files -> "${runDir(meta)}/qc/fastplong" }
    }

    // Assembly //
    assembler_params {
        path { meta, _files -> assemblyDir(meta) }
    }
    contigs {
        // Short-read contigs are renamed to <id>_cleaned.contigs.fa.gz, long-read contigs keep their name
        path { meta, contigs ->
            contigs >> "${assemblyDir(meta)}/${meta.short_reads ? "${meta.id}_cleaned.contigs.fa.gz" : contigs.name}"
        }
    }
    decontamination {
        path { meta, _files -> "${assemblyDir(meta)}/decontamination" }
    }
    coverage {
        path { meta, _files -> "${assemblyDir(meta)}/coverage" }
    }
    quast {
        path { meta, _files -> "${assemblyDir(meta)}/qc/quast" }
    }

    // Reports //
    multiqc_study {
        path { meta, _files -> "${studyDir(meta)}/multiqc" }
    }
    multiqc_run {
        path { meta, _files -> "${assemblyDir(meta)}/qc/multiqc" }
    }
    versions {
        path 'pipeline_info'
    }
    assembled_runs {
        path '.'
    }
    qc_failed_runs {
        path '.'
    }
}
