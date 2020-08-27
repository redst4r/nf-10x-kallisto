# nf-10x-kallisto

This nextflow pipeline can be used to requantify 10x experiments (processed with
cellranger) with kallisto.

## mkfastq_and_kallisto
### Usage
```
nextflow run mkfastq_and_kallisto.nf \
   --bamfile <cellranger_bamfile>
   --kallisto_index <...> \
   --kallisto_gene_map <...> \
   --chemistry <10xv2|10xv3> \
   --barcode_whitelist <...> \
   --outdir <...>
```

### Steps:
1. `mkfastq` to turn the bam-file into fastqs
2. `kallisto bus` to pseudo-align the reads
3. `bustools correct` to correct the cell-barcodes and sort them
4. `bustools count` for gene-wise and equivalence-wise counting of expression


### Output
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
    ├── bustools_metrics
    │   └── bus_output.json
    └── sort_bus
        └── bus_output
            ├── matrix.ec
            ├── output.corrected.sort.bus
            ├── run_info.json
            └── transcripts.txt
```

Note that `kallisto/sort_bus/output.corrected.sort.bus` will be pretty big (a few GB)

## kallisto_from_fastq.nf 
### Usage
Similar to mkfastq_and_kallisto, but starting from fastq files already

```
nextflow run fastq_and_kallisto.nf \
   --readsglob <something_L00*_R{1,2}.fastq.gz> \
   --outdir <directory> \
   --chemistry <10xv2|10xv3> \
```

### Output
Same as mkfastq_and_kallisto.nf


## mkfastq.nf
Just turns a 10x/cellranger bamfile into the fastqs. 

### Usage
```
nextflow run mkfastq.nf \
  --bamfile <.bam>\
   --outdir <directory> \
  --publish_mode <copy|symlink>
```
`--publish_mode copy`  creates copies of the files in outdir, `symlink` just links from the nextflow tmp-directory (to save space)
