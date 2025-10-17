library(here)
library(Seurat)
library(SpatialExperiment)
library(scater)
library(dplyr)
library(tidyverse)
library(harmony)
library(scry)
library(schex)
library(BayesSpace)

setClassUnion("ExpData", c("matrix", "SummarizedExperiment"))

data_name <- "spatialLIBD_DLPFC_12_3_7_12"
load(file = here(data_name, "results", paste0(data_name, "_spe_qc.Rdata")))

svgs_before <- read.csv(file = here(data_name, "results", paste0(data_name, "_svgs.csv")),row.names=1)

set.seed(230)
spe <- scater::runPCA(spe, ncomponents = 50,
                      subset_row = svgs_before$gene_id,
                      BSPARAM = BiocSingular::RandomParam())

spe <- RunHarmony(spe, "subject")
spe <- runUMAP(spe, dimred = "HARMONY", name = "UMAP-HARMONY")

colData(spe)$row <- spe$array_row
colData(spe)$col <- spe$array_col

metadata(spe)$BayesSpace.data <- list(platform = "Visium", is.enhanced = FALSE)

k=7
message("Running spatialCluster()")
Sys.time()
spe <- spatialCluster(spe, use.dimred = "HARMONY", q = k,nrep=10000)
Sys.time()

bayesSpace_name <- paste0("bayesSpace_captureArea_", k)
colnames(colData(spe))[ncol(colData(spe))] <- bayesSpace_name

save(spe, file = here(data_name, "results", paste0(data_name, "_spe_harmony_BayesSpace_preBatchSVG_7.Rdata")))
