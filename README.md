# nf-10x-kallisto

This nextflow pipeline can be used to requantify 10x experiments (processed with
cellranger) with kallisto.

## Usage

```
nextflow run mkfastq_and_kallisto.nf \
   --bamfile <cellranger_bamfile>
   --kallisto_index <...> \
   --kallisto_gene_map <...> \
   --chemistry <10xv2|10xv3> \
   --barcode_whitelist <...> \
   --outdir <...>
```

## Steps:
1. `mkfastq` to turn the bam-file into fastqs
2. `kallisto bus` to pseudo-align the reads
3. `bustools correct` to correct the cell-barcodes and sort them
4. `bustools count` for gene-wise and equivalence-wise counting of expression


## Output
Here's what will be contained in the output folder:

```
.
└── kallisto
    ├── bustools_counts
    │   ├── bus_output_eqcount
    │   │   ├── tcc.barcodes.txt
    │   │   ├── tcc.ec.txt
    │   │   └── tcc.mtx
    │   └── bus_output_genecount
    │       ├── gene.barcodes.txt
    │       ├── gene.genes.txt
    │       └── gene.mtx
    └── bustools_metrics
        └── bus_output.json
```
