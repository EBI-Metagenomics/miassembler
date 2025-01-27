## Introduction

The **ebi-metagenomics/miassembler** supports long reads assembly with [Flye](https://github.com/mikolmogorov/Flye). 

Reads are initially quality-checked with [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/) before and after a [fastp](https://github.com/OpenGene/fastp) evaluation check. The initial QC step also involves decontamination from the human genome and a possible accessory host or environmental genome with [Minimap2](https://github.com/lh3/minimap2). 

The output from fastp is used to determine whether data are low-quality or high-quality, i.e. if Q20 is greater or smaller than 0.8. Based on this value and the sequencing platform used on the original sample, a different branch of the miassembler will be triggered and will execute the following steps:

1. Oxford Nanopore Technology + Low Quality data: 
    - **Adapters trimming** with [Canu](https://canu.readthedocs.io/en/latest/)
    - **Assembly** with Flye `--nano-raw`
    - **Assembly polish** with [racon](https://github.com/isovic/racon)
    - **Assembly polish** with [medaka](https://github.com/nanoporetech/medaka)

2. Oxford Nanopore Technology + High Quality data: 

    -  **Adapters trimming** with [Porechop_ABI](https://github.com/bonsai-team/Porechop_ABI)
    - **Assembly** with Flye `--nano-hq`

3. PacBio + Low Quality data: 

    -  **Adapters trimming** with Canu
    - **Assembly** with Flye `--pacbio-raw`
    - **Assembly polish** with racon

4. PacBio + High Quality data (HiFi): 

    -  **Adapters trimming** with [HiFiAdapterFilt](https://github.com/sheinasim-USDA/HiFiAdapterFilt)
    - **Assembly** with Flye `--pacbio-hifi`

All branches then converge into a post-assembly decontamination step where contigs are newly checked against the human genome and the potential host. 

For low quality ONT and PacBio data, [proovframe](https://github.com/thackl/proovframe) is applied to contigs to apply frameshift correction to contigs. 

Finally, assembly coverage is calculated using [MetaBAT2](https://bitbucket.org/berkeleylab/metabat/src/master/) metabat2_jgi_summarizebamcontigdepths for per contig depth and [Samtools idxstats](http://www.htslib.org/doc/samtools-idxstats.html) for alignment summary statistics.


### Required DBs:

- `--reference_genome`: reference genome in FASTA format
- `--diamond_db`: pre-formatted diamond db to use in the frameshift correction step. Internally, we use NCBI_nr. 