library(SpatialExperiment)
library(here)
library(nnSVG)
library(dplyr)
library(scater)
library(scran)

setClassUnion("ExpData", c("matrix", "SummarizedExperiment"))

data_name <- "spatialLIBD_DLPFC_12_3_7_12"
load(file = here(data_name, "results", paste0(data_name, "_spe_qc.Rdata")))

# run nnSVG filtering for mitochondrial genes and low-expressed genes
spe <- filter_genes(
  spe,
  filter_genes_ncounts = 2,
  filter_genes_pcspots = 0.6,
  filter_mito=F
)

colData(spe)$sample_id <- factor(colData(spe)$sample_id)

sample_ids <- levels(colData(spe)$sample_id)

res_list <- as.list(rep(NA, length(sample_ids)))
names(res_list) <- sample_ids

for (s in seq_along(sample_ids)) {

    # select sample
    ix <- colData(spe)$sample_id == sample_ids[s]
    spe_sub <- spe[, ix]

    dim(spe_sub)

    # remove any zeros introduced by filtering
    ix_zeros <- colSums(counts(spe_sub)) == 0
    if (sum(ix_zeros) > 0) {
        spe_sub <- spe_sub[, !ix_zeros]
    }

    print(dim(spe_sub))

    # re-calculate logcounts after filtering
    spe_sub <- computeLibraryFactors(spe_sub)
    spe_sub <- logNormCounts(spe_sub)

    # run nnSVG
    set.seed(123)
    spe_sub <- nnSVG(spe_sub,n_threads=15,verbose=F)

    # store results for this sample
    res_list[[s]] <- rowData(spe_sub)
    rm(spe_sub)
}

save(res_list, file=here(data_name,"results","res_list_intermediate_expanded.Rdata"))

# match results from each sample and store in matching rows
res_ps <- matrix(NA, nrow = nrow(spe), ncol = length(sample_ids))
rownames(res_ps) <- rownames(spe)
colnames(res_ps) <- sample_ids


for (s in seq_along(sample_ids)) {
  stopifnot(colnames(res_ps)[s] == sample_ids[s])
  stopifnot(colnames(res_ps)[s] == names(res_list)[s])
  
  rownames_s <- rownames(res_list[[s]])
  res_ps[rownames_s, s] <- res_list[[s]][, "padj"]
}

n_signif <- apply(res_ps, 1, function(r) sum(r < 0.05, na.rm = TRUE))


# remove genes that were filtered out in all samples
ix_allna <- apply(res_ps, 1, function(r) all(is.na(r)))
res_ps <- res_ps[!ix_allna, ]

dim(res_ps)

# summary table
df_summary <- data.frame(
    gene_id = rownames(res_ps),
    gene_name = rowData(spe)[rownames(res_ps), "gene_name"],
    n_signif = unname(n_signif),
    row.names = rownames(res_ps)
)

# summary table of "replicated" SVGs (i.e. genes that are highly ranked in at
# least x samples)
df_summaryReplicated <- df_summary[df_summary$n_signif >= 2, ]

dim(df_summaryReplicated)

head(df_summaryReplicated)

# remove mito genes
is_mito <- grepl("(^MT-)|(^mt-)", df_summaryReplicated$gene_name)
df_summaryReplicated <- df_summaryReplicated[!is_mito,]

write.csv(df_summaryReplicated, here(data_name, "results", paste0(data_name, "_expanded_svgs.csv")))
