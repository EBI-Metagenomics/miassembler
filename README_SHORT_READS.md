## Introduction

**ebi-metagenomics/miassembler** is a bioinformatics pipeline for the assembly of short metagenomic reads using [SPAdes](https://doi.org/10.1089/cmb.2012.0021) or [MEGAHIT](https://doi.org/10.1093/bioinformatics/btv033).

1. Read QC using [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)
2. Present QC for raw reads and assembly [MultiQC](http://multiqc.info/)
3. Performs assembly using [MEGAHIT](https://github.com/voutcn/megahit) and [SPAdes](http://cab.spbu.ru/software/spades/), and checks assembly quality using [Quast](http://quast.sourceforge.net/quast)
4. Removes contaminated contigs using [BLASTN](https://blast.ncbi.nlm.nih.gov/Blast.cgi?CMD=Web&PAGE_TYPE=Blastdocs) and [SeqKit](https://bioinf.shenwei.me/seqkit/)
5. Calculates assembly coverage using [MetaBAT2](https://bitbucket.org/berkeleylab/metabat/src/master/) metabat2_jgi_summarizebamcontigdepths for per contig depth and [Samtools idxstats](http://www.htslib.org/doc/samtools-idxstats.html) for alignment summary statistics.

### Required DBs:

- `--blast_reference_genomes_folder`: mandatory **human_phiX** is provided on [FTP](https://ftp.ebi.ac.uk/pub/databases/metagenomics/pipelines/references/)
- `--bwamem2_reference_genomes_folder`: mandatory **human_phiX** is provided on [FTP](https://ftp.ebi.ac.uk/pub/databases/metagenomics/pipelines/references/)

Blast and bwa-mem2 reference databases can be generated for any reference genome to polish input sequences with.

#### BWA-MEM2

As explained in [bwa-mem2's README](https://github.com/bwa-mem2/bwa-mem2?tab=readme-ov-file#getting-started):

```
# Use precompiled binaries (recommended)
curl -L https://github.com/bwa-mem2/bwa-mem2/releases/download/v2.2.1/bwa-mem2-2.2.1_x64-linux.tar.bz2 \
  | tar jxf -

# Index your reference genome with
bwa-mem2-2.2.1_x64-linux/bwa-mem2 index ref.fa
```

This will generate multiple index files in a folder. The folder containing them is the one to use as `bwamem2_reference_genomes_folder`.

#### BLAST

```
makeblastdb -in <ref.fa> -dbtype nucl -out <my_db_file>
```

As with bwa-mem2, numerous files will be generated in the same folder, which should be used for `blast_reference_genomes_folder`.