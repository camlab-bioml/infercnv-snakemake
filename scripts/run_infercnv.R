library(infercnv)
library(SingleCellExperiment)
library(tibble)
library(dplyr)
library(tidyr)
library(argparse)

parser <- ArgumentParser()

parser$add_argument("--sce", type="character")
parser$add_argument("--annotation_column", type="character")
parser$add_argument("--sample_column", type="character")
parser$add_argument("--sample", type="character")
parser$add_argument("--tumour_type", type="character")
parser$add_argument("--normal_type", type="character")
parser$add_argument("--cutoff", type="double")
parser$add_argument("--denoise", type="logical")
parser$add_argument("--genelist", type="character")
parser$add_argument("--out_dir", type="character")
parser$add_argument("--threads", type="integer")
parser$add_argument("--leiden_res", type="double")

args <- parser$parse_args()




set.seed(123L)

sce <- readRDS(args$sce)

sub_sample <- sce[, sce$sample == args$sample]

## Some sanity checks
stopifnot(ncol(sub_sample) > 10)
stopifnot(args$annotation_column %in% names(colData(sce)))
stopifnot(args$sample_column %in% names(colData(sce)))

sce_tumour <- sub_sample[,colData(sub_sample)[[args$annotation_column]] == args$tumour_type]

sce_normal <- sce[,colData(sub_sample)[[args$annotation_column]] == args$normal_type]

if(ncol(sce_normal) > 2000) {
    sce_normal <- sce_normal[, sample(ncol(sce_normal), 2000)]
}

patient_sce <- cbind(sce_tumour, sce_normal)

cts <- assay(patient_sce, 'counts')

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

# cell_types <- colData(patient_sce)[,args$annotation_column, drop=FALSE] |> 
#     as.data.frame() |>
#     rownames_to_column("cell_id")

# rem_cell_types <- group_by(cell_types, cell_type) |>
#     tally() |> 
#     filter(n < 2) |>
#     pull(cell_type)

# rem_cells <- cell_types |>
#     filter(cell_type %in% rem_cell_types) |>
#     pull(cell_id)

# cts <- select(cts, -all_of(rem_cells))
# cell_types <- filter(cell_types, !(cell_id %in% rem_cells)) |>
#     column_to_rownames("cell_id")

normal_cell_types <- args$normal_type #  unique(cell_types$cell_type)
# normal_cell_types <- normal_cell_types[normal_cell_types != "Tumour epithelial"]

infer <- CreateInfercnvObject(cts, 
    args$genelist, 
    c(args$normal_type, args$tumour_type), 
    ref_group_names = args$normal_type)

infercnv_obj <- infercnv::run(infer, 
    cutoff = args$cutoff, 
    out_dir = args$out_dir,
    cluster_by_groups=TRUE, 
    plot_steps=FALSE, 
    denoise = args$denoise 
    # noise_filter=snakemake@params$noise_filter, 
    HMM_report_by = c('subcluster'), 
    tumor_subcluster_pval = 0.01,
    HMM=TRUE, 
    no_prelim_plot=TRUE, 
    num_threads=args$threads,
    leiden_resolution = args$leiden_res, 
    png_res=360)