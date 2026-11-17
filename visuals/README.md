The metro schema for this repo was generated from a simplified version of the dag file generated using the following command:

```bash
nextflow run main.nf -preview -with-dag dag.mmd
```

The dag.mmd was given to Claude alongside the miassembler.nf, short_reads_assembler.nf, and long_reads_assembler.nf files to collapse the subworkflows and decycle the graph. The new simplified dag file was used as input for [nf-metro](https://seqeralabs.github.io/nf-metro/latest/playground/) where labels were fixed in an interactive way. The style brand is Seqera. The final source and figures were exported from nf-metro and uploaded to the GitHub repo at visuals dir.
