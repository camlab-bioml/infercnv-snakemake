library(infercnv)
library(SingleCellExperiment)
library(tibble)
library(dplyr)
library(tidyr)

set.seed(123L)

sce <- readRDS(snakemake@input$sce)

sub_sample <- sce[, sce$sample == snakemake@wildcards$sample]

sce_tumour <- sub_sample[,sub_sample$cell_type == snakemake@params$tumour_type]

sce_normal <- sce[,sub_sample$cell_type == snakemake@params$normal_cell_types]

if(ncol(sce_normal) > 2000) {
    sce_normal <- sce_normal[, sample(ncol(sce_normal), 2000)]
}

patient_sce <- cbind(tumour, normal)

cts <- assays(patient_sce, 'counts')

cts_df <- cts |> 
    as.matrix() |>
    as.data.frame() |> 
    rownames_to_column("gene") |>
    mutate(gene = gsub("ENSG[0-9]*-", "", gene))

cts <- cts_df |>
    pivot_longer(-gene, names_to = "cell_id", values_to = "counts") |>
    group_by(gene, cell_id) |>
    summarize(counts = sum(counts)) |>
    pivot_wider(names_from = "cell_id", values_from = "counts") |>
    as.data.frame() |>
    column_to_rownames("gene")

cell_types <- colData(patient_sce)[,snakemake@params$annotation_column, drop=FALSE] |> 
    as.data.frame() |>
    rownames_to_column("cell_id")

rem_cell_types <- group_by(cell_types, cell_type) |>
    tally() |> 
    filter(n < 2) |>
    pull(cell_type)

rem_cells <- cell_types |>
    filter(cell_type %in% rem_cell_types) |>
    pull(cell_id)

cts <- select(cts, -all_of(rem_cells))
cell_types <- filter(cell_types, !(cell_id %in% rem_cells)) |>
    column_to_rownames("cell_id")

normal_cell_types <- unique(cell_types$cell_type)
normal_cell_types <- normal_cell_types[normal_cell_types != "Tumour epithelial"]

infer <- CreateInfercnvObject(cts, 
    snakemake@input$genelist, 
    cell_types, 
    ref_group_names = snakemake@params$normal_type)

infercnv_obj <- infercnv::run(infer, 
    cutoff = 0, 
    out_dir = snakemake@params$out_dir, 
    cluster_by_groups=TRUE, 
    plot_steps=FALSE, 
    denoise = snakemake@params$denoise, 
    noise_filter=snakemake@params$noise_filter, 
    HMM_report_by = c('subcluster'), 
    tumor_subcluster_pval = 0.01,
    HMM=TRUE, 
    no_prelim_plot=TRUE, 
    num_threads=snakemake@threads, 
    png_res=360)