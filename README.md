## Introduction

**ebi-metagenomics/miassembler** is a bioinformatics pipeline for the assembly of short metagenomic reads using [SPAdes](https://doi.org/10.1089/cmb.2012.0021) or [MEGAHIT](https://doi.org/10.1093/bioinformatics/btv033).

1. Read QC using [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)
2. Present QC for raw reads and assembly [MultiQC](http://multiqc.info/)
3. Performs assembly using [MEGAHIT](https://github.com/voutcn/megahit) and [SPAdes](http://cab.spbu.ru/software/spades/), and checks assembly quality using [Quast](http://quast.sourceforge.net/quast)
4. Removes contaminated contigs using [BLASTN](https://blast.ncbi.nlm.nih.gov/Blast.cgi?CMD=Web&PAGE_TYPE=Blastdocs) and [SeqKit](https://bioinf.shenwei.me/seqkit/)
5. Calculates assembly coverage using [MetaBAT2](https://bitbucket.org/berkeleylab/metabat/src/master/) metabat2_jgi_summarizebamcontigdepths for per contig depth and [Samtools idxstats](http://www.htslib.org/doc/samtools-idxstats.html) for alignment summary statistics.

This pipeline is still in early development. It's mostly a direct port of the mi-automation assembly generation pipeline. Some of the bespoke scripts used to remove contaminated contigs or to calculate the coverage of the assembly were replaced with tools provided by the community ([SeqKit](https://doi.org/10.1371/journal.pone.0163962) and [quast](https://doi.org/10.1093/bioinformatics/btu153) respectively).

> [!NOTE]
> This pipeline uses the [nf-core](https://nf-co.re) template with some tweaks, but it's not part of nf-core.

## Usage

> [!WARNING]
> It only runs in Codon using Slurm ATM.

Pipeline help:

```bash
Typical pipeline command:

  nextflow run ebi-metagenomics/miassembler --help

Input/output options
  --study_accession                       [string]  The ENA Study secondary accession
  --reads_accession                       [string]  The ENA Run primary accession
  --private_study                         [boolean] To use if the ENA study is private
  --samplesheet                           [string]  Path to comma-separated file containing information about the raw reads with the prefix to be used.
  --assembler                             [string]  The short reads assembler (accepted: spades, metaspades, megahit)
  --single_end                            [boolean] Force the single_end value for the study / reads
  --library_strategy                      [string]  Force the library_strategy value for the study / reads (accepted: metagenomic, metatranscriptomic,
                                                    genomic, transcriptomic, other)
  --library_layout                        [string]  Force the library_layout value for the study / reads (accepted: single, paired)
  --spades_version                        [string]  null [default: 3.15.5]
  --megahit_version                       [string]  null [default: 1.2.9]
  --reference_genome                      [string]  The genome to be used to clean the assembly, the genome will be taken from the Microbiome Informatics
                                                    internal directory (accepted: chicken.fna, salmon.fna, cod.fna, pig.fna, cow.fna, mouse.fna,
                                                    honeybee.fna, rainbow_trout.fna, ...)
  --blast_reference_genomes_folder        [string]  The folder with the reference genome blast indexes, defaults to the Microbiome Informatics internal
                                                    directory.
  --bwamem2_reference_genomes_folder      [string]  The folder with the reference genome bwa-mem2 indexes, defaults to the Microbiome Informatics internal
                                                    directory.
  --remove_human_phix                     [boolean] Remove human and phiX reads pre assembly, and contigs matching those genomes. [default: true]
  --human_phix_blast_index_name           [string]  Combined Human and phiX BLAST db. [default: human_phix]
  --human_phix_bwamem2_index_name         [string]  Combined Human and phiX bwa-mem2 index. [default: human_phix]
  --min_contig_length                     [integer] Minimum contig length filter. [default: 500]
  --min_contig_length_metatranscriptomics [integer] Minimum contig length filter for metaT. [default: 200]
  --assembly_memory                       [integer] Default memory allocated for the assembly process. [default: 100]
  --spades_only_assembler                 [boolean] Run SPAdes/metaSPAdes without the error correction step. [default: true]
  --outdir                                [string]  The output directory where the results will be saved. You have to use absolute paths to storage on Cloud
                                                    infrastructure. [default: results]
  --email                                 [string]  Email address for completion summary.
  --multiqc_title                         [string]  MultiQC report title. Printed as page header, used for filename if not otherwise specified.

Generic options
  --multiqc_methods_description           [string]  Custom MultiQC yaml file containing HTML including a methods description.
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

### Samplesheet

The samplesheet is a comma-separated file (.csv) with the following columns:

- study_accession: Unique identifier for the study.
- reads_accession: Unique identifier for the reads.
- fastq_1: Full path to the first FastQ file.
- fastq_2: Full path to the second FastQ file (for paired-end reads). Leave empty if single-end.
- library_layout: Either single or paired.
- library_strategy: One of metagenomic, metatranscriptomic, genomic, transcriptomic, or other.

The header row is mandatory. Below is an example of a minimum required samplesheet:

```csv
study_accession,reads_accession,fastq_1,fastq_2,library_layout,library_strategy
PRJ1,ERR1,/path/to/reads/ERR1_1.fq.gz,/path/to/reads/ERR1_2.fq.gz,paired,metagenomic
PRJ2,ERR2,/path/to/reads/ERR2.fq.gz,,single,genomic
```

#### Full Samplesheet Example

The pipeline can handle both single-end and paired-end reads. A full samplesheet for different library layouts and strategies might look like this:

```csv
study_accession,reads_accession,fastq_1,fastq_2,library_layout,library_strategy,assembler
PRJ1,ERR1,/path/to/reads/ERR1_1.fq.gz,/path/to/reads/ERR1_2.fq.gz,paired,metagenomic
PRJ2,ERR2,/path/to/reads/ERR2.fq.gz,,single,genomic,metaspades
PRJ3,ERR3,/path/to/reads/ERR3_1.fq.gz,/path/to/reads/ERR3_2.fq.gz,paired,transcriptomic
```

It also allow to specify which assembler and how much memory in GB to use for the assembly process. Assembler and assembler memory columns:

- assembler: One of spades, metaspades, or megahit.
- assembly_memory: Integer value specifying the memory allocated for the assembly process.

Example with additional columns:

```csv
study_accession,reads_accession,fastq_1,fastq_2,library_layout,library_strategy,assembler,assembly_memory
PRJ1,ERR1,/path/to/reads/ERR1_1.fq.gz,/path/to/reads/ERR1_2.fq.gz,paired,metagenomic,spades,16
PRJ2,ERR2,/path/to/reads/ERR2.fq.gz,,single,genomic,megahit,32
```

## Outputs

The outputs of the pipeline are organized as follows:

```
results
├── pipeline_info
├── DRP0076
│   └── DRP007622
│       ├── DRR2807
│       │   └── DRR280712
│       │       ├── assembly
│       │       │   └── megahit
│       │       │       └── 1.2.9
│       │       │           ├── coverage
│       │       │           ├── decontamination
│       │       │           └── qc
│       │       │               ├── multiqc
│       │       │               └── quast
│       │       │                   └── DRR280712
│       │       └── qc
│       │           ├── fastp
│       │           └── fastqc
│       └── multiqc
└── SRP1154
    └── SRP115494
        ├── multiqc
        ├── SRR5949
        │   └── SRR5949318
        │       ├── assembly
        │       │   └── metaspades
        │       │       └── 3.15.5
        │       │           ├── coverage
        │       │           ├── decontamination
        │       │           └── qc
        │       │               ├── multiqc
        │       │               └── quast
        │       │                   └── SRR5949318
        │       └── qc
        │           ├── fastp
        │           └── fastqc
        └── SRR6180
            └── SRR6180434 --> QC Failed (not assembled)
                └── qc
                    ├── fastp
                    └── fastqc
```

The nested structure based on ENA Study and Reads accessions was created to suit the Microbiome Informatics team’s needs. The benefit of this structure is that results from different runs of the same study won’t overwrite any results.

### Top Level Reports

#### MultiQC

The pipeline produces two [MultiQC](https://multiqc.info) reports: one per study and one per run. These reports aggregate statistics related to raw reads, read QC, assembly, and assembly QC.

The run-level MultiQC report is generated for runs that passed QC and were assembled. The study-level MultiQC report includes all runs; however, runs without assemblies will not have assembly stats included.

#### QC failed runs

QC failed runs are filtered out to prevent downstream assembly failures.

Runs that fail QC checks are excluded from the assembly process. These runs are listed in the file `qc_failed_runs.csv`, along with the corresponding exclusion message. Assembling such runs may cause the pipeline to fail or produce very poor assemblies.

Example:

```csv
SRR6180434,filter_ratio_threshold_exceeded
```

##### Runs exclusion messages

| Exclusion Message                 | Description                                                                                                                                                                                                                                                                            |
| --------------------------------- | -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `filter_ratio_threshold_exceeded` | The maximum allowed ratio of reads after filtering compared to before filtering. If exceeded, it flags excessive filtering. The default value is 0.9, meaning that if more than 90% of the reads are filtered out, the threshold is considered exceeded, and the run is not assembled. |
| `low_reads_count_threshold`       | The minimum number of reads required after filtering. If below, it flags a low read count, and the run is not assembled.                                                                                                                                                               |

#### Assembled Runs

Runs that were successfully assembled are listed in a CSV file named `assembled_runs.csv`. This file contains the run accession, assembler, and assembler version used.

Example:

```csv
DRR280712,megahit,1.2.9
SRR5949318,metaspades,3.15.5
```

## Tests

There is a very small test data set ready to use:

```bash
nextflow run main.nf -resume -profile test,docker
```

It's also possible to run the [nf-test](https://www.nf-test.com/) suite with

```bash
nf-test test
```

### End to end tests

Two end-to-end tests can be launched (with megahit and metaspades) with the following command:

```bash
pytest tests/workflows/ --verbose
```
