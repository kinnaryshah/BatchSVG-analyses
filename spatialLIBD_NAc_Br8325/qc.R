library(here)
library(SpatialExperiment)
library(nnSVG)
library(scran)
library(DropletUtils)
library(scuttle)
library(spatialLIBD)
library(tidyverse)
library(jaffelab)
library(sessioninfo)

setClassUnion("ExpData", c("matrix", "SummarizedExperiment"))

data_name <- "spatialLIBD_NAc_Br8325"

load(file=here(data_name,"results",paste0(data_name,"_spe_all.Rdata")))

spe <- spe[,colData(spe)$sample_id %in% c("V13F06-313_C1", "V13M06-376_A1", "V13M06-376_B1", "V13M06-376_C1", "V13M06-376_D1")]

print(dim(spe))

# keep only spots over tissue
spe <- spe[,colData(spe)[,"in_tissue"]==1]

# remove genes with zero expression
ix_zero_genes <- rowSums(counts(spe)) == 0
table(ix_zero_genes)
if (sum(ix_zero_genes) > 0) {
    spe <- spe[!ix_zero_genes, ]
}

# remove spots with zero expression
ix_zero_spots <- colSums(counts(spe)) == 0
if (sum(ix_zero_spots) > 0) {
    spe <- spe[, !ix_zero_spots]
}

print(dim(spe))

save(spe, file=here("data_name","results", paste0(data_name,"_spe_qc.Rdata")))
