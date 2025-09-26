library(SpatialExperiment)
library(spatialLIBD)
library(BatchSVG)
library(dplyr)
library(ggplot2)
library(scran)
library(ComplexHeatmap)
library(scuttle)
library(scry)
library(tidyverse)

featureSelect <- function(input, batch_effects = NULL, VGs = NULL,
                          verbose = TRUE) {

  if (!is.null(batch_effects)) {
    stopifnot(is.character(batch_effects))
    for (batch in batch_effects) {
      if (!batch %in% names(colData(input))) {
        stop("The batch_effect is not a valid column")
      }
    }
  } else {
    stop("Please provide a valid batch_effect.")
  }
  
  input <- input[rowData(input)$gene_id %in% VGs, ]
  rownames(input) <- rowData(input)$gene_id
  
  if (verbose) message("Running feature selection with batch...")
  bd <- devianceFeatureSelection(input, assay = "counts", fam = "binomial")
  bd_df <- cbind.data.frame("gene_id"=rownames(bd),
                            "gene_name"=rowData(bd)$gene_name, "dev"= rowData(bd)$binomial_deviance,
                            "rank"=(nrow(bd) + 1) - rank(rowData(bd)$binomial_deviance))
  rownames(bd_df) <- bd_df$gene
  
  results_list <- vector("list", length(batch_effects))
  names(results_list) <- batch_effects  
  
  for (i in seq_along(batch_effects)) {
    batch <- batch_effects[i]
    
    if (verbose) message("Batch Effect: ", batch)
    if (verbose) message("Running feature selection without batch...")
    
    batch_data <- colData(input)[[batch]]
    bd_batch <- devianceFeatureSelection(input, assay = "counts",
                                         fam = "binomial",batch = as.factor(batch_data))
    
    bd_batch_df <- cbind.data.frame(
      "gene_id"=rownames(bd_batch),
      "gene_name"=rowData(bd_batch)$gene_name,
      "dev"= rowData(bd_batch)$binomial_deviance,
      "rank"=(nrow(bd_batch)+1)-rank(rowData(bd_batch)$binomial_deviance))
    rownames(bd_batch_df) <- bd_batch_df$gene
    
    if (verbose) message("Calculating deviance and rank difference...")
    
    batch_df <- left_join(bd_df, bd_batch_df, by=c("gene_id", "gene_name"),
                          suffix=c("_default",paste0("_", batch)))
    
    batch_df$d_diff <- 
      (batch_df$dev_default- batch_df[[paste0("dev_", batch)]]) /
      batch_df[[paste0("dev_", batch)]]
    batch_df[[paste0("nSD_dev_", batch)]] <- 
      (batch_df$d_diff - mean(batch_df$d_diff)) / sd(batch_df$d_diff)
    
    batch_df$r_diff <- 
      batch_df[[paste0("rank_", batch)]] - batch_df$rank_default
    batch_df[[paste0("nSD_rank_", batch)]] <- 
      (batch_df$r_diff - mean(batch_df$r_diff)) / sd(batch_df$r_diff)
    
    results_list[[i]] <- batch_df
  }
  results_list
}


set.seed(123)

sce_sn = fetch_data(type="spatialDLPFC_snRNAseq")
unzip(sce_sn, exdir = tempdir())
sce_sn <- HDF5Array::loadHDF5SummarizedExperiment(
  file.path(tempdir(), "sce_DLPFC_annotated")
)
sce_sn
colData(sce_sn)

sce <- sce_sn

# just choose 2 M and 2 F samples
sce <- sce[, which(sce$Sample %in% c("Br2720_mid","Br2743_mid","Br8325_mid","Br6522_mid"))]

# remove genes with zero expression
ix_zero_genes <- rowSums(counts(sce)) == 0
table(ix_zero_genes)
if (sum(ix_zero_genes) > 0) {
  sce <- sce[!ix_zero_genes, ]
}

# remove nuclei with zero expression
ix_zero_nuclei <- colSums(counts(sce)) == 0
if (sum(ix_zero_nuclei) > 0) {
  sce <- sce[, !ix_zero_nuclei]
}

print(dim(sce))

# get HVGs
dec <- modelGeneVar(sce)
hvgs_idx <- getTopHVGs(dec, prop=0.1, row.names = F)
hvgs <- rowData(sce)$gene_id[hvgs_idx]
str(hvgs)

# run BatchSVG
list_batch_df <- featureSelect(input = sce, 
                               batch_effect = c("Sample"), VGs = hvgs)
head(list_batch_df[["Sample"]])

plots <- svg_nSD(list_batch_df = list_batch_df, 
                 sd_interval_dev = 6, sd_interval_rank = 6)
plots

bias_both <- biasDetect(list_batch_df = list_batch_df, threshold = "rank",
           nSD_dev = 6, nSD_rank = 6, plot_point_shape = 23, plot_palette = "RdPu",
           plot_text_size = 4)
bias_both[["Sample"]][["Plot"]]
genes <- bias_both[["Sample"]][["Table"]]$gene_name


bias_both <- biasDetect(list_batch_df = list_batch_df, threshold = "dev",
                        nSD_dev = 7, nSD_rank = 7, plot_point_shape = 23, plot_palette = "RdPu",
                        plot_text_size = 4)
bias_both[["Sample"]][["Plot"]]
genes <- append(genes,bias_both[["Sample"]][["Table"]]$gene_name)
genes <- unique(genes)

# vis in heatmap
log_mat <- logcounts(sce)[genes, , drop = FALSE]

cell_types <- colData(sce)$Sample

long_expr <- as.data.frame(as.matrix(log_mat)) %>%
  rownames_to_column("gene") %>%
  tidyr::pivot_longer(-gene, names_to = "cell", values_to = "logcounts") %>%
  mutate(cell_type = cell_types[match(cell, colnames(sce))])

mean_expr <- long_expr %>%
  group_by(gene, cell_type) %>%
  summarise(mean_logcounts = mean(logcounts), .groups = "drop") %>%
  tidyr::pivot_wider(names_from = cell_type, values_from = mean_logcounts)

heatmap_mat <- as.matrix(column_to_rownames(mean_expr, var = "gene"))

Heatmap(
  heatmap_mat,
  cluster_rows = F, cluster_columns = F, row_names_side = "left",
  show_row_names = TRUE, show_column_names = TRUE,
  heatmap_legend_param = list(
    title = "mean\nlogcount\nexpression", at = c(0, 4),
    labels = c("0", "4")
  ),
)

