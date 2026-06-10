setwd("Desktop/Research/BatchSVG-analyses/")
library(here)
library(SpatialExperiment)
library(scater)
library(ggplot2)
library(ggspavis)
library(spatialLIBD)
library(PCAtools)
library(dplyr)
library(cowplot)
library(bluster)
library(escheR)
library(patchwork)
library(RColorBrewer)
library(tibble)
library(ComplexHeatmap)
library(glmpca)
library(scry)
library(schex)
library(mclust)
library(tidyr)
library(aricode)
library(BatchSVG)

load(file=here("humanHippocampus2024","results", "humanHippocampus2024_sample_id_feat_sel.Rdata"))
vec <- list_batch_df$sample_id$nSD_dev_sample_id[list_batch_df$sample_id$nSD_dev_sample_id > 0]
x_values <- 1:33

counts <- sapply(x_values, function(x) sum(vec > x))

df <- data.frame(
  threshold = x_values,
  count = counts
)

p1 <- ggplot(df, aes(x = threshold, y = count)) +
  geom_point(color = "blue", size = 2) +
  labs(x = "Threshold Value",
       y = "Number of SVGs > Threshold",
       title = "SVGs Exceeding Threshold - HPC Relative Change in Deviance") +
  theme_minimal()

vec <- list_batch_df$sample_id$nSD_rank_sample_id[list_batch_df$sample_id$nSD_rank_sample_id > 0]
x_values <- 1:17

counts <- sapply(x_values, function(x) sum(vec > x))

df <- data.frame(
  threshold = x_values,
  count = counts
)

p2 <- ggplot(df, aes(x = threshold, y = count)) +
  geom_point(color = "blue", size = 2) +
  labs(x = "Threshold Value",
       y = "Number of SVGs > Threshold",
       title = "SVGs Exceeding Threshold - HPC Rank Difference") +
  theme_minimal()


load(file=here("spatialLIBD_DLPFC_12_3_7_12_expanded","results", "spatialLIBD_DLPFC_12_3_7_12_expanded_subject_feat_sel.Rdata"))
vec <- list_batch_df$subject$nSD_dev_subject[list_batch_df$subject$nSD_dev_subject > 0]
x_values <- 1:32

counts <- sapply(x_values, function(x) sum(vec > x))

df <- data.frame(
  threshold = x_values,
  count = counts
)

p3 <- ggplot(df, aes(x = threshold, y = count)) +
  geom_point(color = "blue", size = 2) +
  labs(x = "Threshold Value",
       y = "Number of SVGs > Threshold",
       title = "SVGs Exceeding Threshold - dlPFC Relative Change in Deviance") +
  theme_minimal()

vec <- list_batch_df$subject$nSD_rank_subject[list_batch_df$subject$nSD_rank_subject > 0]
x_values <- 1:33

counts <- sapply(x_values, function(x) sum(vec > x))

df <- data.frame(
  threshold = x_values,
  count = counts
)

p4 <- ggplot(df, aes(x = threshold, y = count)) +
  geom_point(color = "blue", size = 2) +
  labs(x = "Threshold Value",
       y = "Number of SVGs > Threshold",
       title = "SVGs Exceeding Threshold - dlPFC Rank Difference") +
  theme_minimal()

data_name = "humanHippocampus2024"
png(here(data_name,"plots","sensitivity.png"),height=8,width=12,unit="in",res=300)

wrap_plots(p1,p2,p3,p4,nrow=2)

dev.off()
