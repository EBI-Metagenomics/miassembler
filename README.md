## Introduction

**ebi-metagenomics/miassembler** is a bioinformatics pipeline for the assembly of long and short metagenomic reads.

This pipeline supports both short and long reads; however, it does not yet support hybrid assemblies.

The steps of the pipeline for short- and long-reads processing are outlined in the following documents:

    [Short Reads](README_SHORT_READS.md)
    [Long Reads](README_LONG_READS.md)

This pipeline is mostly a direct port of the mi-automation assembly generation pipeline. Some of the bespoke scripts used to remove contaminated contigs or to calculate the coverage of the assembly were replaced with tools provided by the community ([SeqKit](https://doi.org/10.1371/journal.pone.0163962) and [quast](https://doi.org/10.1093/bioinformatics/btu153) respectively).

> [!NOTE]
> This pipeline uses the [nf-core](https://nf-co.re) template with some tweaks, but it's not part of nf-core.

## Usage

Pipeline help:

```angular2html
Typical pipeline command:

  nextflow run main.nf -profile <docker/singularity/...> --samplesheet samplesheet.csv --outdir <OUTDIR>

--help                                    [boolean, string] Show the help message for all top level parameters. When a parameter is given to `--help`, the full help message of that parameter will be printed.

Input/output options
  --samplesheet                           [string]  Path to comma-separated file containing information about the raw reads with the prefix (read accession) to be used.
  --study_accession                       [string]  The ENA Study secondary accession
  --reads_accession                       [string]  The ENA Run primary accession
  --private_study                         [boolean] To use if the ENA study is private, *this feature ony works on EBI infrastructure at the moment*
  --reference_genome                      [string]  The genome to be used to clean the assembly, the genome will be taken from the Microbiome Informatics internal directory  (accepted: chicken.fna, salmon.fna, cod.fna, pig.fna,
cow.fna, mouse.fna, honeybee.fna, rainbow_trout.fna, rat.fna, sheep.fna, soybean.fna, zebrafish.fna)
  --reference_genomes_folder              [string]  The folder with the reference genomes, defaults to the Microbiome Informatics internal directory.
  --blast_reference_genomes_folder        [string]  The folder with the reference genome blast indexes, defaults to the Microbiome Informatics internal directory.
  --bwamem2_reference_genomes_folder      [string]  The folder with the reference genome bwa-mem2 indexes, defaults to the Microbiome Informatics internal directory.
  --diamond_db                            [string]  Path to diamond db (e.g. NCBI-nr) to perform frameshift correction [default: ]
  --outdir                                [string]  The output directory where the results will be saved. You have to use absolute paths to storage on Cloud infrastructure. [default: results]
  --email                                 [string]  Email address for completion summary.
  --multiqc_title                         [string]  MultiQC report title. Printed as page header, used for filename if not otherwise specified.

Input metadata options
  --assembler                             [string]  The short or long reads assembler  (accepted: spades, metaspades, megahit, flye)
  --long_reads_assembler_config           [string]  Configuration to use flye with.  (accepted: nano-raw, nano-corr, nano-hq, pacbio-raw, pacbio-corr, pacbio-hifi)
  --flye_version                          [string]  [default: 2.9]
  --spades_version                        [string]  [default: 3.15.5]
  --megahit_version                       [string]  [default: 1.2.9]
  --remove_human_phix                     [boolean] Set true to removing human and phiX reads pre-assembly, and contigs matching those genomes. [default: true]
  --remove_human                          [boolean] Set true to removing human reads pre-assembly, and contigs matching those genomes. [default: true]
  --human_phix_blast_index_name           [string]  Filename of the combined Human and phiX BLAST db. [default: human_phix]
  --human_phix_bwamem2_index_name         [string]  Filename of the combined Human and phiX bwa-mem2 index. [default: human_phix]
  --human_fasta_prefix                    [string]  Prefix in the filename of the human genome reference. [default: human]
  --assembly_memory                       [number]  Default memory allocated for the assembly process. [default: 100]
  --spades_only_assembler                 [boolean] Run SPAdes/metaSPAdes without the error correction step. [default: true]

Reads QC options
  --short_reads_filter_ratio_threshold    [number]  The maximum fraction of reads that are allowed to be filtered out. If exceeded, it flags excessive filtering. The default value is 0.1, meaning that if less than 10% of the reads
are retained after filtering, the threshold is considered exceeded, and the run is not assembled. [default: 0.1]
  --short_reads_low_reads_count_threshold [number]  The minimum number of reads required after filtering. If below, it flags a low read count and the run is not assembled. [default: 1000]
  --long_reads_min_read_length            [integer] Minimum read length for pre-assembly quality filtering [default: 200]
  --long_reads_pacbio_quality_threshold   [number]  The Q20 threshold that a pacbio sample needs to exceed to be labelled as HiFi. [default: 0.8]

Assembly QC options
  --short_reads_min_contig_length         [integer] Minimum contig length filter for short reads. [default: 500]
  --short_reads_min_contig_length_metat   [integer] Minimum contig length filter for short reads metaT. [default: 200]
  --short_reads_contig_threshold          [integer] Minimum number of contigs in human+phiX+host cleaned assembly. [default: 2]

Generic options
  --multiqc_methods_description           [string] Custom MultiQC yaml file containing HTML including a methods description.
```

You can run this pipeline with two options:

### Command-line parameters

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

```bash
nextflow run ebi-metagenomics/miassembler \
  -profile codon_slurm \
  --samplesheet tests/samplesheet/test.csv
```

The samplesheet is a comma-separated file (.csv) with the following columns:

- study_accession: Unique identifier for the study.
- reads_accession: Unique identifier for the reads.
- fastq_1: Full path to the first FastQ file.
- fastq_2: Full path to the second FastQ file (for paired-end reads). Leave empty if single-end.
- library_layout: Either single or paired.
- library_strategy: One of metagenomic, metatranscriptomic, genomic, transcriptomic, or other.
- platform: Relevant for long reads, requiring either ont or pb for nanopore or pacbio, respectively.
- assembler: Relevant for short reads, where either megahit, metaspades, or spades can be picked. Flye is also supported
- assembly_memory: Integer value specifying the memory allocated for the assembly process.
- assembler_config: Configuration to use flye with. One of "nano-raw", "nano-corr", "nano-hq", "pacbio-raw", "pacbio-corr", "pacbio-hifi".

The header row is mandatory.

#### Full Samplesheet Example

The pipeline can handle both single-end and paired-end reads. A full samplesheet for different library layouts and strategies might look like this:

```csv
study_accession,reads_accession,fastq_1,fastq_2,library_layout,library_strategy,assembler
PRJ1,ERR1,/path/to/reads/ERR1_1.fq.gz,/path/to/reads/ERR1_2.fq.gz,paired,metagenomic
PRJ2,ERR2,/path/to/reads/ERR2.fq.gz,,single,genomic,metaspades
PRJ3,ERR3,/path/to/reads/ERR3_1.fq.gz,/path/to/reads/ERR3_2.fq.gz,paired,transcriptomic
```

Example with additional columns:

```csv
study_accession,reads_accession,fastq_1,fastq_2,library_layout,library_strategy,assembler,assembly_memory,assembler_config,platform
PRJ1,ERR1,/path/to/reads/ERR1_1.fq.gz,/path/to/reads/ERR1_2.fq.gz,paired,metagenomic,spades,16
PRJ2,ERR2,/path/to/reads/ERR2.fq.gz,,single,genomic,flye,32,nano-hq,ont
```

### ENA Private Data

The pipeline includes a module to download private data from ENA using the EMBL-EBI FIRE (File Replication) system. This system is restricted for use within the EMBL-EBI network and will not work unless connected to that network.

If you have private data to assemble, you must provide the full path to the files on a system that Nextflow can access.

#### Microbiome Informatics Team

To process private data, the pipeline should be launched with the `--private_study` flag, and the samplesheet must include the private FTP (transfer services) paths. The `download_from_fire` module will be utilized to download the files.

This module uses [Nextflow secrets](https://www.nextflow.io/docs/latest/secrets.html#how-it-works). Specifically, it requires the `FIRE_ACCESS_KEY` and `FIRE_SECRET_KEY` secrets to authenticate and download the files.

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

### Coverage

The pipeline reports the coverage values for the assembly using two mechanisms: `jgi_summarize_bam_contig_depths` and a custom whole assembly coverage and coverage depth.

#### jgi_summarize_bam_contig_depths

This tool summarizes the depth of coverage for each contig from BAM files containing the mapped reads. It quantifies the extent to which contigs in an assembly are covered by these reads. The output is a tabular file, with rows representing contigs and columns displaying the summarized coverage values from the BAM files. This summary is useful for binning contigs or estimating abundance in various metagenomic datasets.

This file is generated per assembly and stored in the following location (e.g., for study `SRP115494` and run `SRR6180434`): `SRP1154/SRP115494/multiqc/SRR5949/SRR5949318/assembly/metaspades/3.15.5/coverage/SRR6180434_coverage_depth_summary.tsv.gz`

##### Example output of `jgi_summarize_bam_contig_depths`

| contigName                       | contigLen | totalAvgDepth | SRR6180434_sorted.bam | SRR6180434_sorted.bam-var |
| -------------------------------- | --------- | ------------- | --------------------- | ------------------------- |
| NODE_1_length_539_cov_105.072314 | 539       | 273.694       | 273.694               | 74284.7                   |

###### Explanation of the Columns:

1. **contigName**: The name or identifier of the contig (e.g., `NODE_1_length_539_cov_105.072314`). This is usually derived from the assembly process and may include information such as the contig length and coverage.

2. **contigLen**: The length of the contig in base pairs (e.g., `539`).

3. **totalAvgDepth**: The average depth of coverage across the entire contig from all BAM files (e.g., `273.694`). This represents the total sequencing coverage averaged across the length of the contig. This value will be the same as the sample avg. depth in assemblies of a single sample.

4. **SRR6180434_sorted.bam**: The average depth of coverage for the specific sample represented by this BAM file (e.g., `273.694`). This shows how well the contig is covered by reads.

5. **SRR6180434_sorted.bam-var**: The variance in the depth of coverage for the same BAM file (e.g., `74284.7`). This gives a measure of how uniform or uneven the read coverage is across the contig.

#### Coverage JSON

The pipeline calculates two key metrics: coverage and coverage depth for the entire assembly. The coverage is determined by dividing the number of assembled base pairs by the total number of base pairs before filtering. Coverage depth is calculated by dividing the number of assembled base pairs by the total length of the assembly, provided the assembly length is greater than zero. These metrics provide insights into how well the reads cover the assembly and the average depth of coverage across the assembled contigs. The script that calculates this number is [calculate_assembly_coverage.py](bin/calculate_assembly_coverage.py).

The pipeline creates a JSON file with the following content:

```json
{
  "coverage": 0.04760503915318373,
  "coverage_depth": 273.694
}
```

The file is stored in (e.g. for study `SRP115494` and run `SRR6180434`) -> `SRP1154/SRP115494/multiqc/SRR5949/SRR5949318/assembly/metaspades/3.15.5/coverage/SRR6180434_coverage.json`

### Top Level Reports

#### MultiQC

The pipeline produces two [MultiQC](https://multiqc.info) reports: one per study and one per run. These reports aggregate statistics related to raw reads, read QC, assembly, and assembly QC.

The run-level MultiQC report is generated for runs that passed QC and were assembled. The study-level MultiQC report includes all runs; however, runs without assemblies will not have assembly stats included.

#### QC failed runs

QC failed runs are filtered out to prevent downstream assembly failures.

Runs that fail QC checks are excluded from the assembly process. These runs are listed in the file `qc_failed_runs.csv`, along with the corresponding exclusion message. Assembling such runs may cause the pipeline to fail or produce very poor assemblies.

Example:

```csv
SRR6180434,short_reads_filter_ratio_threshold_exceeded
```

##### Runs exclusion messages

| Exclusion Message                             | Description                                                                                                                                                                                                                                                                          |
| --------------------------------------------- | ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------ |
| `short_reads_filter_ratio_threshold_exceeded` | The maximum fraction of reads that are allowed to be filtered out. If exceeded, it flags excessive filtering. The default value is 0.1, meaning that if less than 10% of the reads are retained after filtering, the threshold is considered exceeded, and the run is not assembled. |
| `short_reads_low_reads_count_threshold`       | The minimum number of reads required after filtering. If below, it flags a low read count, and the run is not assembled.                                                                                                                                                             |
| `short_reads_contig_threshold`                | The minimum number of contigs allowed after human+phiX+host cleaning. If below it flags a low contig count and the cleaned assembly isn't generated.                                                                                                                                 |

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
