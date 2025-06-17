# ebi-metagenomics/miassembler: Documentation

**ebi-metagenomics/miassembler** is a bioinformatics pipeline for the assembly of long and short metagenomic reads.

This pipeline supports both short and long reads; however, it does not yet support hybrid assemblies.

The pipeline steps for both processes follow.

## Short reads

For the assembly of short metagenomic reads the pipeline uses [SPAdes](https://doi.org/10.1089/cmb.2012.0021) or [MEGAHIT](https://doi.org/10.1093/bioinformatics/btv033).

1. Read QC using [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)
2. Present QC for raw reads and assembly [MultiQC](http://multiqc.info/)
3. Performs assembly using [MEGAHIT](https://github.com/voutcn/megahit) and [SPAdes](http://cab.spbu.ru/software/spades/), and checks assembly quality using [Quast](http://quast.sourceforge.net/quast)
4. Removes contaminated contigs using [Minimap2](https://github.com/lh3/minimap2) and [SeqKit](https://bioinf.shenwei.me/seqkit/)
5. Calculates assembly coverage using [MetaBAT2](https://bitbucket.org/berkeleylab/metabat/src/master/) metabat2_jgi_summarizebamcontigdepths for per contig depth and [Samtools idxstats](http://www.htslib.org/doc/samtools-idxstats.html) for alignment summary statistics.

### Required DBs:

- `--reference_genomes_folder`: Path to the folder that contains all required reference genomes in FASTA format as well as their indexes generated with bwa-mem2. **phiX** is provided on [FTP](https://ftp.ebi.ac.uk/pub/databases/metagenomics/pipelines/references/)

bwa-mem2 index files must be generated for all genomes provided as references for read decontamination.

#### How to index reference genome with BWA-MEM2

As explained in [bwa-mem2's README](https://github.com/bwa-mem2/bwa-mem2?tab=readme-ov-file#getting-started):

```
# Use precompiled binaries (recommended)
curl -L https://github.com/bwa-mem2/bwa-mem2/releases/download/v2.2.1/bwa-mem2-2.2.1_x64-linux.tar.bz2 \
  | tar jxf -

# Index your reference genome with
bwa-mem2-2.2.1_x64-linux/bwa-mem2 index ref.fa
```

This will generate multiple index files in a folder. Place them together with source FASTA file to `reference_genomes_folder`.

## Long reads

The pipeline uses [Flye](https://github.com/mikolmogorov/Flye) to assemble long reads.

Reads are initially quality-checked with [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/) before and after a [fastp](https://github.com/OpenGene/fastp) evaluation check. The initial QC step also involves decontamination from the human genome and a possible accessory host or environmental genome with [Minimap2](https://github.com/lh3/minimap2).

The output from fastp is used to determine whether data are low-quality or high-quality, i.e. if Q20 is greater or smaller than 0.8. Based on this value and the sequencing platform used on the original sample, a different branch of the miassembler will be triggered and will execute the following steps:

1. Oxford Nanopore Technology + Low Quality data:

   - **Adapters trimming** with [Canu](https://canu.readthedocs.io/en/latest/)
   - **Assembly** with Flye `--nano-raw`
   - **Assembly polish** with [racon](https://github.com/isovic/racon)
   - **Assembly polish** with [medaka](https://github.com/nanoporetech/medaka)

2. Oxford Nanopore Technology + High Quality data:

   - **Adapters trimming** with [Porechop_ABI](https://github.com/bonsai-team/Porechop_ABI)
   - **Assembly** with Flye `--nano-hq`

3. PacBio + Low Quality data:

   - **Adapters trimming** with Canu
   - **Assembly** with Flye `--pacbio-raw`
   - **Assembly polish** with racon

4. PacBio + High Quality data (HiFi):

   - **Adapters trimming** with [HiFiAdapterFilt](https://github.com/sheinasim-USDA/HiFiAdapterFilt)
   - **Assembly** with Flye `--pacbio-hifi`

All branches then converge into a post-assembly decontamination step where contigs are newly checked against the human genome and the potential host.

For low quality ONT and PacBio data, [proovframe](https://github.com/thackl/proovframe) is applied to contigs to apply frameshift correction to contigs.

Finally, assembly coverage is calculated using [MetaBAT2](https://bitbucket.org/berkeleylab/metabat/src/master/) metabat2_jgi_summarizebamcontigdepths for per contig depth and [Samtools idxstats](http://www.htslib.org/doc/samtools-idxstats.html) for alignment summary statistics.

### Required DBs:

- `--reference_genomes_folder`: Path to the folder that contains all required reference genomes in FASTA format.
- `--diamond_db`: pre-formatted diamond db to use in the frameshift correction step. Internally, we use NCBI_nr.
