## Introduction

**ebi-metagenomics/miassembler** is a bioinformatics pipeline for the assembly of short metagenomic reads using [SPAdes](https://doi.org/10.1089/cmb.2012.0021) or [MEGAHIT](https://doi.org/10.1093/bioinformatics/btv033).

1. Read QC using [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)
2. Present QC for raw reads and assembly [MultiQC](http://multiqc.info/)
3. Performs assembly using [MEGAHIT](https://github.com/voutcn/megahit) and [SPAdes](http://cab.spbu.ru/software/spades/), and checks assembly quality using [Quast](http://quast.sourceforge.net/quast)
4. Removes contaminated contigs using [BLASTN](https://blast.ncbi.nlm.nih.gov/Blast.cgi?CMD=Web&PAGE_TYPE=Blastdocs) and [SeqKit](https://bioinf.shenwei.me/seqkit/)
5. Calculates assembly coverage using [MetaBAT2](https://bitbucket.org/berkeleylab/metabat/src/master/) metabat2_jgi_summarizebamcontigdepths for per contig depth and [Samtools idxstats](http://www.htslib.org/doc/samtools-idxstats.html) for alignment summary statistics.

This pipeline is still in early development. It's mostly a direct port of the mi-automation assembly generation pipeline. Some of the bespoke scripts used to remove contaminated contigs or to calculate the coverage of the assembly were replaced with tools provided by the community ([SeqKit](https://doi.org/10.1371/journal.pone.0163962) and [quast](https://doi.org/10.1093/bioinformatics/btu153) respectively).

## Usage

> [!WARNING]
> It only runs in Codon using Slurm ATM.

Pipeline help:

```bash
nextflow run ebi-metagenomics/miassembler --help

Input/output options
  --study_accession                  [string]  The ENA Study secondary accession
  --reads_accession                  [string]  The ENA Run primary accession
  --private_study                    [boolean] To use if the ENA study is private [default: false]
  --assembler                        [string]  The short reads assembler (accepted: spades, metaspades, megahit) [default: metaspades]
  --reference_genome                 [string]  The genome to be used to clean the assembly, the genome will be taken from the Microbiome Informatics internal
                                               directory (accepted: chicken.fna, salmon.fna, cod.fna, pig.fna, cow.fna, mouse.fna, honeybee.fna,
                                               rainbow_trout.fna, rat.fna, ...)
  --blast_reference_genomes_folder   [string]  The folder with the reference genome blast indexes, defaults to the Microbiome Informatics internal
                                               directory.
  --bwamem2_reference_genomes_folder [string]  The folder with the reference genome bwa-mem2 indexes, defaults to the Microbiome Informatics internal
                                               directory.
  --remove_human_phix                [boolean] Remove human and phiX reads pre assembly, and contigs matching those genomes. [default: true]
  --human_phix_blast_index_name      [string]  Combined Human and phiX BLAST db. [default: human_phix]
  --human_phix_bwamem2_index_name    [string]  Combined Human and phiX bwa-mem2 index. [default: human_phix]
  --min_contig_length                [integer] Minimum contig length filter. [default: 500]
  --assembly_memory                  [integer] Default memory allocated for the assembly process. [default: 100]
  --spades_only_assembler            [boolean] Run SPAdes/metaSPAdes without the error correction step. [default: true]
  --outdir                           [string]  The output directory where the results will be saved. You have to use absolute paths to storage on Cloud
                                               infrastructure.
  --email                            [string]  Email address for completion summary.
  --multiqc_title                    [string]  MultiQC report title. Printed as page header, used for filename if not otherwise specified.

Generic options
  --multiqc_methods_description      [string]  Custom MultiQC yaml file containing HTML including a methods description.
```

Example:

```bash
nextflow run ebi-metagenomics/miassembler \
  -profile codon_slurm \
  --assembler metaspades \
  --reference_genome human \
  --outdir testing_results \
  --study_accession SRP002480 \
  --reads_accession SRR1631361
```

## Outputs

The outputs of the pipeline are organized as follows:

```
results/SRP1154
└── SRP115494
    └── SRR6180
        └── SRR6180434
            ├── assembly
            │   └── metaspades
            │       └── 3.15.5
            │           ├── coverage
            │           ├── decontamination
            │           └── qc
            │               ├── multiqc
            │               └── quast
            └── qc
                ├── fastp
                └── fastqc

```

The nested structure based on ENA Study and Reads accessions was created to suit the Microbiome Informatics team’s needs. The benefit of this structure is that results from different runs of the same study won’t overwrite any results.

## Tests

There is a very small test data set ready to use:

```bash
nextflow run main.nf -resume -profile test,docker
```

### End to end tests

Two end-to-end tests can be launched (with megahit and metaspades) with the following command:

```bash
pytest tests/workflows/ --verbose
```
