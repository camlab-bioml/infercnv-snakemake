
## Step 1: install correct snakemake version

```bash
pipenv install
pipenv shell
```

## Step 2: edit config file

```
input_rds: peng-scRNASeq-manually-filtered-sce-tumour-normal-assigned.rds
annotation_column: cell_type
sample_column: donor

tumour_type: Tumour epithelial
normal_type: Normal epithelial

cutoff: 0.1
denoise: TRUE

samples:
  - CRR034499
  - CRR034500
```

Key entries:

- `input_rds`

## Step 3: Run snakemake
