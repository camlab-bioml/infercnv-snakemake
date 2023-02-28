
## Step 1: Install relevant Python packages

```bash
pipenv install --python 3.8
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

- `input_rds`: input SingleCellExperiment containing both tumour and normal. Gene names should either be `[SYMBOL]` or `[ENSEMBLID]-[SYMBOL]`
- `annotation_column`: what column of `colData(sce)` says which cells are tumour vs normal?
- `sample_column`: what column of `colData(sce)` refers to sample (patient/donor)? The pipeline runs inferCNV once per donor
- `tumour_type` The ID in `colData(sce)[[annotation_column]]` that refers to the cells that are tumour/malignant
- `normal_type` The ID in `colData(sce)[[annotation_column]]` that refers to the cells that are 'normal'
- `cutoff`,`denoise` parameters passed to inferCNV
- `samples` list of samples (patients/donors/etc) to run inferCNV over. Should be present in `colData(sce)[[sample_column]]`

## Step 3: Run the pipeline

```bash
snakemake -j1 --configfile config/peng-test.yml --use-singularity --singularity-args "--bind /home/campbell/share/:/home/campbell/share/"
```

replacing `/home/campbell/share` with your directory.

## Step 4: Inspect results


# Notes

- This pipeline aggregates genes to gene symbol by summing