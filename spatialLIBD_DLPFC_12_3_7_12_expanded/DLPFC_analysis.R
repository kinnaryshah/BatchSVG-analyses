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

load(file=here("spatialLIBD_DLPFC_12_3_7_12_expanded","results", "spatialLIBD_DLPFC_12_3_7_12_expanded_spe_PRECAST_preBatchSVG_6.Rdata"))
spe_pre <- spe
spe_pre$PRECAST_cluster <- unfactor(spe_pre$PRECAST_cluster)
spe_pre$PRECAST_cluster[spe_pre$PRECAST_cluster == 1] <- "L1"
spe_pre$PRECAST_cluster[spe_pre$PRECAST_cluster == 2] <- "L3/4/5/6 (1)"
spe_pre$PRECAST_cluster[spe_pre$PRECAST_cluster == 3] <- "L3/4/5/6 (2)"
spe_pre$PRECAST_cluster[spe_pre$PRECAST_cluster == 4] <- "L2"
spe_pre$PRECAST_cluster[spe_pre$PRECAST_cluster == 5] <- "WM (2)"
spe_pre$PRECAST_cluster[spe_pre$PRECAST_cluster == 6] <- "WM"

load(file=here("spatialLIBD_DLPFC_12_3_7_12_expanded","results", "spatialLIBD_DLPFC_12_3_7_12_expanded_spe_PRECAST_postBatchSVG_6.Rdata"))
spe_post <- spe
spe_post$PRECAST_cluster <- unfactor(spe_post$PRECAST_cluster)
spe_post$PRECAST_cluster[spe_post$PRECAST_cluster == 1] <- "WM"
spe_post$PRECAST_cluster[spe_post$PRECAST_cluster == 2] <- "L5/6"
spe_post$PRECAST_cluster[spe_post$PRECAST_cluster == 3] <- "L1"
spe_post$PRECAST_cluster[spe_post$PRECAST_cluster == 4] <- "L2"
spe_post$PRECAST_cluster[spe_post$PRECAST_cluster == 5] <- "L3/4"
spe_post$PRECAST_cluster[spe_post$PRECAST_cluster == 6] <- "WM (2)"

svgs <- read.csv(here("spatialLIBD_DLPFC_12_3_7_12_expanded","results","spatialLIBD_DLPFC_12_3_7_12_expanded_svgs.csv"), row.names=1,
                 check.names=F)

filt_svgs <- read.csv(here("spatialLIBD_DLPFC_12_3_7_12_expanded","results","spatialLIBD_DLPFC_12_3_7_12_expanded_subject_7_7_filt_svgs.csv"),
                      check.names=F)

# marker gene heatmaps
source("../dlpfc_genes.r")
dlpfc.genes$Inhb <- NULL

rownames(spe_pre) <- rowData(spe_pre)$gene_name

genes <- intersect(unlist(dlpfc.genes),rownames(spe_pre))
cell_type_per_gene <- c(rep("Micro.Vasc",4),rep("Astro",4),
                        rep("L2",5),rep("L3",3),rep("L4",4),
                        rep("L5",4),rep("L6",3),rep("Oligo",2))

spe_pre$PRECAST_cluster <- as.factor(spe_pre$PRECAST_cluster)

log_mat <- logcounts(spe_pre)[genes, , drop = FALSE]
cell_types <- colData(spe_pre)$PRECAST_cluster

long_expr <- as.data.frame(as.matrix(log_mat)) %>%
  rownames_to_column("gene") %>%
  tidyr::pivot_longer(-gene, names_to = "cell", values_to = "logcounts") %>%
  mutate(cell_type = cell_types[match(cell, colnames(spe_pre))])

mean_expr <- long_expr %>%
  group_by(gene, cell_type) %>%
  summarise(mean_logcounts = mean(logcounts), .groups = "drop") %>%
  tidyr::pivot_wider(names_from = cell_type, values_from = mean_logcounts)

heatmap_mat <- as.matrix(column_to_rownames(mean_expr, var = "gene"))

# reorder based on genes variable
heatmap_mat <- heatmap_mat[genes, ]

heatmap_mat <- t(scale(t(heatmap_mat)))

colors <- brewer.pal(n = 8, name = "Dark2")
names(colors) <- unique(cell_type_per_gene)

data_name = "spatialLIBD_DLPFC_12_3_7_12_expanded"
png(here(data_name,"plots","PRECAST_cluster_pre_marker_heatmap_logcounts_DLPFC.png"),height=5,width=8,unit="in",res=300)

Heatmap(
  t(heatmap_mat),
  cluster_rows = F, cluster_columns = F, row_names_side = "left",
  show_row_names = TRUE, show_column_names = TRUE,
  heatmap_legend_param = list(
    title = "scaled avg.\nlogcount\nexpression", at = c(-4, 0, 4),
    labels = c("-4", "0", "4")
  ),
  top_annotation = HeatmapAnnotation(
    " " = cell_type_per_gene,
    col = list(" " = colors),
    show_legend = T
  ),
)

dev.off()


rownames(spe_post) <- rowData(spe_post)$gene_name

genes <- intersect(unlist(dlpfc.genes),rownames(spe_post))
cell_type_per_gene <- c(rep("Micro.Vasc",4),rep("Astro",4),
                        rep("L2",5),rep("L3",3),rep("L4",4),
                        rep("L5",4),rep("L6",3),rep("Oligo",2))

spe_post$PRECAST_cluster <- as.factor(spe_post$PRECAST_cluster)

log_mat <- logcounts(spe_post)[genes, , drop = FALSE]
cell_types <- colData(spe_post)$PRECAST_cluster

long_expr <- as.data.frame(as.matrix(log_mat)) %>%
  rownames_to_column("gene") %>%
  tidyr::pivot_longer(-gene, names_to = "cell", values_to = "logcounts") %>%
  mutate(cell_type = cell_types[match(cell, colnames(spe_post))])

mean_expr <- long_expr %>%
  group_by(gene, cell_type) %>%
  summarise(mean_logcounts = mean(logcounts), .groups = "drop") %>%
  tidyr::pivot_wider(names_from = cell_type, values_from = mean_logcounts)

heatmap_mat <- as.matrix(column_to_rownames(mean_expr, var = "gene"))

# reorder based on genes variable
heatmap_mat <- heatmap_mat[genes, ]

heatmap_mat <- t(scale(t(heatmap_mat)))

colors <- brewer.pal(n = 8, name = "Dark2")
names(colors) <- unique(cell_type_per_gene)

data_name = "spatialLIBD_DLPFC_12_3_7_12_expanded"
png(here(data_name,"plots","PRECAST_cluster_post_marker_heatmap_logcounts_DLPFC.png"),height=5,width=8,unit="in",res=300)

Heatmap(
  t(heatmap_mat),
  cluster_rows = F, cluster_columns = F, row_names_side = "left",
  show_row_names = TRUE, show_column_names = TRUE,
  heatmap_legend_param = list(
    title = "scaled avg.\nlogcount\nexpression", at = c(-4, 0, 4),
    labels = c("-4", "0", "4")
  ),
  top_annotation = HeatmapAnnotation(
    " " = cell_type_per_gene,
    col = list(" " = colors),
    show_legend = T
  ),
)

dev.off()


# set colors for each set of clusters
colors_pre <- brewer.pal(n = 6, name = "Paired")
colors_pre <- setNames(colors_pre, c("WM","WM (2)","L2","L3/4/5/6 (2)","L3/4/5/6 (1)","L1"))

colors_post <- brewer.pal(n = 6, name = "Paired")
colors_post <- setNames(colors_post, c("WM","WM (2)","L2","L5/6","L3/4","L1"))

# vis PRECAST clusters

data_name = "spatialLIBD_DLPFC_12_3_7_12_expanded"
png(here(data_name,"plots","PRECAST_cluster_7_pre_clusters.png"),height=5,width=6,unit="in",res=300)

p <- plotCoords(spe_pre,sample_id="sample_id",annotate = "PRECAST_cluster",assay_name = "logcounts",
                pal=colors_pre) +
  ggtitle("Domains With All SVGs") +
  labs(color = "PRECAST Cluster")

p +
  facet_wrap(
    "sample_id",
    nrow = 1
  )


dev.off()

data_name = "spatialLIBD_DLPFC_12_3_7_12_expanded"
png(here(data_name,"plots","PRECAST_cluster_7_post_clusters.png"),height=5,width=6,unit="in",res=300)

p <- plotCoords(spe_post,sample_id="sample_id",annotate = "PRECAST_cluster",assay_name = "logcounts",
                pal=colors_post) +
  ggtitle("Domains Without BatchSVGs") +
  labs(color = "PRECAST Cluster")


p +
  facet_wrap(
    "sample_id",
    nrow = 1
  )
dev.off()





genes_exp <- c("")
clusters_list_pre <- c()

for (gene in genes_exp) {
  print(gene)
  
  spe_pre_1 <- spe_pre[, which(spe_pre$sample_id == "151509")]
  p1 <- make_escheR(spe_pre_1) |> add_fill(var="PRECAST_cluster", point_size = 1) |> add_ground(var="PRECAST_cluster", stroke=0.1, point_size = 1) +
    scale_fill_manual(values=colors_pre) +
    scale_color_manual(values=colors_pre) +
    theme(legend.position="none") +
    ggtitle("Sample 151509") +
    labs(fill="Before BatchSVG") + guides(color="none")
  
  spe_pre_2 <- spe_pre[, which(spe_pre$sample_id == "151671")]
  p2 <- make_escheR(spe_pre_2) |> add_fill(var="PRECAST_cluster", point_size = 1) |> add_ground(var="PRECAST_cluster", stroke=0.1, point_size = 1) +
    scale_fill_manual(values=colors_pre) +
    scale_color_manual(values=colors_pre) +
    theme(legend.position="none") +
    ggtitle("Sample 151671") +
    labs(fill="Before BatchSVG") + guides(color="none")
  
  spe_pre_3 <- spe_pre[, which(spe_pre$sample_id == "151676")]
  p3 <- make_escheR(spe_pre_3) |> add_fill(var="PRECAST_cluster", point_size = 1) |> add_ground(var="PRECAST_cluster", stroke=0.1, point_size = 1) +
    scale_color_manual(values=colors_pre) +
    scale_fill_manual(values=colors_pre) +
    ggtitle("Sample 151676") +
    labs(fill="Before BatchSVG") + guides(color="none")
  
  clusters_list_pre <- c(clusters_list_pre, list(p1,p2,p3))
  
}


genes_exp <- c("")
clusters_list_post <- c()

for (gene in genes_exp) {
  print(gene)
  
  spe_post_1 <- spe_post[, which(spe_post$sample_id == "151509")]
  p1 <- make_escheR(spe_post_1) |> add_fill(var="PRECAST_cluster", point_size = 1) |> add_ground(var="PRECAST_cluster", stroke=0.1, point_size = 1) +
    scale_fill_manual(values=colors_post) +
    scale_color_manual(values=colors_post) +
    theme(legend.position="none") +
    labs(fill="After BatchSVG") + guides(color="none")
  
  spe_post_2 <- spe_post[, which(spe_post$sample_id == "151671")]
  p2 <- make_escheR(spe_post_2) |> add_fill(var="PRECAST_cluster", point_size = 1) |> add_ground(var="PRECAST_cluster", stroke=0.1, point_size = 1) +
    scale_fill_manual(values=colors_post) +
    scale_color_manual(values=colors_post) +
    theme(legend.position="none")  +
    labs(fill="After BatchSVG") + guides(color="none")
  
  spe_post_3 <- spe_post[, which(spe_post$sample_id == "151676")]
  p3 <- make_escheR(spe_post_3) |> add_fill(var="PRECAST_cluster", point_size = 1) |> add_ground(var="PRECAST_cluster", stroke=0.1, point_size = 1) +
    scale_color_manual(values=colors_post) +
    scale_fill_manual(values=colors_post)  +
    labs(fill="After BatchSVG") + guides(color="none")
  
  
  clusters_list_post <- c(clusters_list_post, list(p1,p2,p3))
  
}


data_name = "spatialLIBD_DLPFC_12_3_7_12_expanded"
png(here(data_name,"plots","PRECAST_cluster_both_clusters.png"),height=9,width=9,unit="in",res=300)

wrap_plots(clusters_list_pre[[1]],clusters_list_pre[[2]],clusters_list_pre[[3]],
           clusters_list_post[[1]],clusters_list_post[[2]],clusters_list_post[[3]],
           ncol=3) +
  plot_layout(guides = 'collect', axes = "collect")

dev.off()


# vis QC metrics breakdown by old and new clusters

spe_pre <- scuttle::addPerCellQC(
  spe_pre,
  subsets = list(Mito = grep("^MT-", rowData(spe_pre)$gene_name)),
  BPPARAM = BiocParallel::MulticoreParam(4)
)

spe_post <- scuttle::addPerCellQC(
  spe_post,
  subsets = list(Mito = grep("^MT-", rowData(spe_post)$gene_name)),
  BPPARAM = BiocParallel::MulticoreParam(4)
)

p1 <- plotColData(spe_pre, x = "PRECAST_cluster", y = "subsets_Mito_percent", colour_by = "PRECAST_cluster") +
  scale_color_manual(values=colors_pre) +
  xlab("") +
  ylab("Mito Percent") +
  labs(colour = "PRECAST Domain\nBefore BatchSVG") +
  theme(
    axis.ticks.x=element_blank())

p2 <- plotColData(spe_post, x = "PRECAST_cluster", y = "subsets_Mito_percent", colour_by = "PRECAST_cluster") +
  scale_color_manual(values=colors_post) +
  xlab("") +
  ylab("Mito Percent") +
  labs(colour = "PRECAST Domain\nAfter BatchSVG") +
  theme(
    axis.ticks.x=element_blank())

p3 <- plotColData(spe_pre, x = "PRECAST_cluster", y = "sum", colour_by = "PRECAST_cluster") +
  scale_color_manual(values=colors_pre) +
  xlab("") +
  ylab("Sum") +
  labs(colour = "PRECAST Domain\nBefore BatchSVG") +
  theme(
    axis.ticks.x=element_blank())

p4 <- plotColData(spe_post, x = "PRECAST_cluster", y = "sum", colour_by = "PRECAST_cluster") +
  scale_color_manual(values=colors_post) +
  xlab("") +
  ylab("Sum") +
  labs(colour = "PRECAST Domain\nAfter BatchSVG") +
  theme(
    axis.ticks.x=element_blank())

p5 <- plotColData(spe_pre, x = "PRECAST_cluster", y = "detected", colour_by = "PRECAST_cluster") +
  scale_color_manual(values=colors_pre) +
  xlab("") +
  ylab("Detected") +
  labs(colour = "PRECAST Domain\nBefore BatchSVG") +
  theme(
    axis.ticks.x=element_blank())

p6 <- plotColData(spe_post, x = "PRECAST_cluster", y = "detected", colour_by = "PRECAST_cluster") +
  scale_color_manual(values=colors_post) +
  xlab("") +
  ylab("Detected") +
  labs(colour = "PRECAST Domain\nAfter BatchSVG") +
  theme(
    axis.ticks.x=element_blank())

data_name = "spatialLIBD_DLPFC_12_3_7_12_expanded"

png(here(data_name,"plots","PRECAST_cluster_6_qc_violins.png"),height=5,width=15,unit="in",res=300)

wrap_plots(p1,p3,p5,p2,p4,p6,
           guides="collect") 
dev.off()

# remove NAs to use aricode() package
na_idx <- is.na(spe_pre$cell_type)
NMI(spe_pre$cell_type[!na_idx],spe_pre$cluster[!na_idx])

na_idx <- is.na(spe_post$cell_type)
NMI(spe_post$cell_type[!na_idx],spe_post$cluster[!na_idx])



# sil plots for pre 
# filter only SVGs
spe_pre <- spe_pre[rownames(spe_pre) %in% svgs$gene_id,]

set.seed(9)
message("running nullResiduals - ", Sys.time())
res_pre <- nullResiduals(spe_pre,
                         fam = "poisson",
                         type = "pearson",
                         assay = "counts"   #, batch = res$brain giving errors
)

set.seed(10)
message("running PCA - ", Sys.time())

res_pre <- scater::runPCA(res_pre, ncomponents = 50,
                          exprs_values='poisson_pearson_residuals',
                          scale = TRUE, name = "pp-GLM-PCA",
                          BSPARAM = BiocSingular::RandomParam())

hex <- make_hexbin(res_pre, nbins = 100, dimension_reduction = "pp-GLM-PCA", use_dims = c(1, 2))
label_df <- make_hexbin_label(hex, col = "PRECAST_cluster")

# get % var from 
#reducedDims(res_pre)[["pp-GLM-PCA"]]

data_name = "spatialLIBD_DLPFC_12_3_7_12_expanded"
png(here(data_name,"plots","PRECAST_cluster_pre_PCA.png"),height=9,width=11,unit="in",res=300)

pca1 <- plot_hexbin_meta(hex, col = "PRECAST_cluster", action = "majority", xlab = "PC1 (35.36%)", ylab = "PC2 (14.04%)") +
  ggtitle("Before BatchSVG") + theme(legend.position = "right") +
  scale_fill_manual(values=colors_pre)

pca1

dev.off()

sil.approx <- approxSilhouette(reducedDim(res_pre, "pp-GLM-PCA"), clusters=res_pre$PRECAST_cluster)
sil.approx

sil.data <- as.data.frame(sil.approx)
sil.data$closest <- factor(ifelse(sil.data$width > 0, res_pre$PRECAST_cluster, sil.data$other))
sil.data$cluster <- res_pre$PRECAST_cluster


sil.means <- sil.data %>%
  group_by(cluster) %>%
  summarise(mean_width = mean(width, na.rm = TRUE))

p1 <- ggplot(sil.data, aes(x=cluster, y=width, colour=closest)) +
  ggbeeswarm::geom_quasirandom(method="smiley", alpha=0.5, size=0.3) +
  scale_colour_manual(values=colors_pre) +
  geom_text(
    data = sil.means,
    aes(x = cluster, y = -0.6, label = sprintf("Mean = %.2f", mean_width)),
    inherit.aes = FALSE,
    vjust = 1.2,
    size = 2.8
  ) +
  theme_bw() +
  ylim(-0.6, 0.5) +
  xlab("Cluster") +
  ylab("Silhouette Width (All SVGs)") +
  labs(colour = "Closest Domain\nBefore BatchSVG")


# sil plots for post 

# filter only SVGs
spe_post <- spe_post[rownames(spe_post) %in% filt_svgs$gene_id,]

set.seed(9)
message("running nullResiduals - ", Sys.time())
res_post <- nullResiduals(spe_post,
                          fam = "poisson",
                          type = "pearson",
                          assay = "counts"   
)

set.seed(10)
message("running PCA - ", Sys.time())

res_post <- scater::runPCA(res_post, ncomponents = 50,
                           exprs_values='poisson_pearson_residuals',
                           scale = TRUE, name = "pp-GLM-PCA",
                           BSPARAM = BiocSingular::RandomParam())

hex <- make_hexbin(res_post, nbins = 100, dimension_reduction = "pp-GLM-PCA", use_dims = c(1, 2))
label_df <- make_hexbin_label(hex, col = "PRECAST_cluster")

data_name = "spatialLIBD_DLPFC_12_3_7_12_expanded"
png(here(data_name,"plots","PRECAST_cluster_post_PCA.png"),height=9,width=11,unit="in",res=300)

pca2 <- plot_hexbin_meta(hex, col = "PRECAST_cluster", action = "majority", xlab = "PC1 (35.30%)", ylab = "PC2 (14.01%)") +
  ggtitle("After BatchSVG") + theme(legend.position = "right") +
  scale_fill_manual(values=colors_post) 

pca2

dev.off()

data_name = "spatialLIBD_DLPFC_12_3_7_12_expanded"
png(here(data_name,"plots","PRECAST_cluster_both_PCA.png"),height=12,width=9,unit="in",res=300)

wrap_plots(pca1,pca2,ncol=1)

dev.off()

sil.approx <- approxSilhouette(reducedDim(res_post, "pp-GLM-PCA"), clusters=res_post$PRECAST_cluster)
sil.approx

sil.data <- as.data.frame(sil.approx)
sil.data$closest <- factor(ifelse(sil.data$width > 0, res_post$PRECAST_cluster, sil.data$other))
sil.data$cluster <- res_post$PRECAST_cluster

sil.means <- sil.data %>%
  group_by(cluster) %>%
  summarise(mean_width = mean(width, na.rm = TRUE))


p2 <- ggplot(sil.data, aes(x=cluster, y=width, colour=closest)) +
  ggbeeswarm::geom_quasirandom(method="smiley",alpha=0.5,size=0.3) +
  scale_colour_manual(values=colors_post) +
  geom_text(
    data = sil.means,
    aes(x = cluster, y = -0.6, label = sprintf("Mean = %.2f", mean_width)),
    inherit.aes = FALSE,
    vjust = 1.2,
    size = 2.8
  ) +
  theme_bw() +
  ylim(-0.6,0.5) +
  xlab("Cluster") +
  ylab("Silhouette Width (Without BatchSVGs)") +
  labs(colour = "Closest Domain\nAfter BatchSVG")


png(here(data_name,"plots","PRECAST_cluster_both_sil.png"),height=6,width=9,unit="in",res=300)

wrap_plots(p1,p2,ncol=1, axes="collect_x")

dev.off()



# vis batch-biased SVGs

load(file=here("spatialLIBD_DLPFC_12_3_7_12_expanded","results", "spatialLIBD_DLPFC_12_3_7_12_expanded_spe_qc.Rdata"))
print(dim(spe))

splot1 <- plotCoords(spe, annotate="ENSG00000256618", assay="logcounts", 
                     sample_id="sample_id", point_size=.1,
                     pal=c("white","black")) +
  ggtitle("MTRNR2L1") +
  labs(color="logcounts")

splot1 <- splot1 +
  facet_wrap(
    "sample_id",
    nrow = 1
  )


splot2 <- plotCoords(spe, annotate="ENSG00000255823", assay="logcounts", 
                     sample_id="sample_id", point_size=.1,
                     pal=c("white","black")) +
  ggtitle("MTRNR2L8") +
  labs(color="logcounts")

splot2 <- splot2 +
  facet_wrap(
    "sample_id",
    nrow = 1
  )

data_name = "spatialLIBD_DLPFC_12_3_7_12_expanded"
png(here(data_name,"plots","spot_plot_MTRNR2L1.png"),height=5,width=6,unit="in",res=300)

splot1

dev.off()

data_name = "spatialLIBD_DLPFC_12_3_7_12_expanded"
png(here(data_name,"plots","spot_plot_MTRNR2L8.png"),height=5,width=6,unit="in",res=300)

splot2

dev.off()

load(file=here("spatialLIBD_DLPFC_12_3_7_12_expanded","results", "spatialLIBD_DLPFC_12_3_7_12_expanded_subject_feat_sel.RData"))

# needed to reformat the function a bit

.theme_dev_point_plot <- function(plot, point_size, point_shape) {
  plot +
    geom_point(size = point_size, shape = point_shape) + 
    scale_x_log10() + 
    scale_y_log10() +
    geom_abline(aes(slope = 1, intercept = 0), lty = 2) +
    theme_bw() + 
    theme(legend.position = "right", aspect.ratio = 1,
          plot.title = element_text(size = 12), 
          plot.subtitle = element_text(size = 10, face = "italic"),
          legend.title = element_text(size = 8),
          axis.title.x = element_text(size = 10),
          axis.title.y = element_text(size = 10)) +
    labs(x= "Deviance without batch", y="Deviance with batch", 
         color = "nSD Deviance",
         title = "Deviance Comparison")
}

.theme_rank_point_plot <- function(plot, point_size, point_shape) {
  plot +
    geom_point(size = point_size, shape = point_shape) + 
    scale_y_reverse() + 
    geom_abline(aes(slope = -1, intercept = 0), lty = 2) +
    theme_bw() + 
    theme(legend.position = "right", aspect.ratio = 1,
          plot.title = element_text(size = 12), 
          plot.subtitle = element_text(size = 10, face = "italic"),
          legend.title = element_text(size = 8),
          axis.title.x = element_text(size = 10),
          axis.title.y = element_text(size = 10)) +
    labs(x= "Rank without batch", y="Rank with batch", 
         color = "nSD Rank",
         title = "Rank Comparison")
}


num_batches <- 1
nSD_rank <- 10
nSD_dev <- 7
plot_point_size <- 1
plot_point_shape <- 16
plot_text_size <- 3
plot_palette <- "RdPu"


biased_list <- vector("list", length(list_batch_df))
names(biased_list) <- names(list_batch_df)

dev_sd_plot <- NULL
rank_sd_plot <- NULL

for (i in seq_along(list_batch_df)) {
  batch <- names(list_batch_df)[i]
  batch_df <- list_batch_df[[batch]]
  stopifnot(is.data.frame(batch_df))
  
  if (!is.null(nSD_dev)) {
    sd_dev <- nSD_dev[i]
    dev_colname <- paste0("nSD_dev_",batch)
    
    batch_df$nSD_bin_dev <- cut(abs(batch_df[[dev_colname]]), right = FALSE,
                                breaks=seq(0,max(batch_df[[dev_colname]]) + sd_dev, 
                                           by=sd_dev), include.lowest=TRUE)
    
    col_pal_dev <- brewer.pal(length(unique(batch_df[["nSD_bin_dev"]])), 
                              plot_palette[i])
    col_pal_dev[1] <- "grey"
    
    dev_sd_plot <- ggplot(batch_df,
                          aes(x = .data[["dev_default"]], y = .data[[paste0("dev_", batch)]],
                              color = .data[["nSD_bin_dev"]]))
    dev_sd_plot <- .theme_dev_point_plot(dev_sd_plot,
                                         point_size = plot_point_size[i], point_shape = plot_point_shape[i])+
      scale_color_manual(values=col_pal_dev) +
      labs(subtitle = "")+
      geom_text_repel(
        aes(label = ifelse(.data[[dev_colname]] > sd_dev,
                           .data[["gene_name"]], "")), size = plot_text_size[i], 
        max.overlaps = Inf,
        show.legend = F) +
      theme(legend.position="bottom")
    
    batch_df$dev_outlier <- batch_df$nSD_dev >= sd_dev
  }
  
  if (!is.null(nSD_rank)) {
    sd_rank <- nSD_rank[i]
    rank_colname <- paste0("nSD_rank_", batch)
    
    batch_df$nSD_bin_rank <- cut(abs(batch_df[[rank_colname]]), right=FALSE,
                                 breaks=seq(0,max(batch_df[[rank_colname]]) + sd_rank, 
                                            by=sd_rank),include.lowest=TRUE)
    
    col_pal_rank <- brewer.pal(length(unique(batch_df$nSD_bin_rank)), 
                               plot_palette[i])
    col_pal_rank[1] <- "grey"
    
    rank_sd_plot <- ggplot(batch_df, 
                           aes(x = .data[["rank_default"]],y = .data[[paste0("rank_", batch)]],
                               color = .data[["nSD_bin_rank"]]))
    rank_sd_plot <- .theme_rank_point_plot(rank_sd_plot,
                                           point_size = plot_point_size[i], point_shape = plot_point_shape[i])+
      scale_color_manual(values = col_pal_rank) +
      labs(subtitle = "")+
      geom_text_repel(
        aes(label = ifelse(.data[[rank_colname]] > sd_rank,
                           .data[["gene_name"]], "")), size = plot_text_size[i], 
        max.overlaps = Inf,
        show.legend = F) +
      theme(legend.position="bottom")
    
  }   
  
}

data_name = "spatialLIBD_DLPFC_12_3_7_12_expanded"
png(here(data_name,"plots","bias_features_final.png"),height=3,width=18,unit="in",res=300)

wrap_plots(dev_sd_plot,rank_sd_plot)

dev.off()