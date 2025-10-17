library(here)
library(PRECAST)
library(Seurat)
library(SpatialExperiment)
library(purrr)
library(dplyr)
library(tidyverse)

setClassUnion("ExpData", c("matrix", "SummarizedExperiment"))

data_name <- "spatialLIBD_DLPFC_12_3_7_12"
load(file = here(data_name, "results", paste0(data_name, "_spe_qc.Rdata")))

colnames(spe) <- spe$key

seuList <- unique(spe$subject) |>
    set_names(unique(spe$subject)) |>
    map(.f = function(id) {
        tmp_spe <- spe[, spe$subject == id]

        tmp_spe$row <- tmp_spe$array_row
        tmp_spe$col <- tmp_spe$array_col

        # browser()
        CreateSeuratObject(
            counts=as.matrix(counts(tmp_spe)),
            meta.data=data.frame(colData(tmp_spe)),
            project="dACC")
    })

svgs_before <- read.csv(file = here(data_name, "results", paste0(data_name, "_svgs.csv")),row.names=1)

set.seed(1)
preobj <- CreatePRECASTObject(seuList = seuList, customGenelist = svgs_before$gene_id,
                              premin.spots = 0, premin.features=0, postmin.spots=0, postmin.features=0)
preobj@seulist

PRECASTObj <- AddAdjList(preobj, platform = "Visium")
PRECASTObj <- AddParSetting(PRECASTObj, Sigma_equal = FALSE,  maxIter = 20, verbose = TRUE)
PRECASTObj <- PRECAST(PRECASTObj, K=7)
PRECASTObj <- SelectModel(PRECASTObj)
seuInt <- IntegrateSpaData(PRECASTObj, species = "Human")

# Merge with spe object
cluster_df <- seuInt@meta.data |>
    mutate(cluster = factor(cluster)) |>
    rename_with(~ paste0("PRECAST_", .x)) |>
    rownames_to_column(var = "key")

col_data_df <- colData(spe) |>
    data.frame() |>
    left_join(cluster_df, by="key")

rownames(col_data_df) <- colnames(spe)
colData(spe)$PRECAST_cluster <- col_data_df$PRECAST_cluster

save(spe, file = here(data_name, "results", paste0(data_name, "_spe_PRECAST_preBatchSVG_7.Rdata")))
