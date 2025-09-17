#!/usr/bin/env nextflow

/*
 * download data into SpatialExperiment objects
 */
process download_data {

    publishDir 'results', mode: 'copy'

    input:
        val data_name

    script:
    """
    #!/usr/bin/env Rscript
    library(here)
    library(SpatialExperiment)
    library(nnSVG)
    library(scran) 
    library(DropletUtils)
    library(scuttle)

    spe <- read10xVisium(here("SpatialBenchVisium/preprocessed/OCT_CytAssist_WT_709/outs/"),
                     type = "sparse",   # use sparse (not HDF5) format
                     data = "raw",     
                     images = "lowres", # specify which image(s) to include
                     load = TRUE)      # specify whether or not to load image(s)


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
    colData(spe)[,"sample_id"] <- rep("OCT_CytAssist_WT_709",dim(spe)[2])
   
    spe1 <- spe

    spe <- read10xVisium(here("SpatialBenchVisium/preprocessed/FFPE_2_WT_709/outs/"),
                     type = "sparse",   # use sparse (not HDF5) format
                     data = "raw",     
                     images = "lowres", # specify which image(s) to include
                     load = TRUE)      # specify whether or not to load image(s)


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
    colData(spe)[,"sample_id"] <- rep("FFPE_2_WT_709",dim(spe)[2])
   
    spe2 <- spe

    spe <- read10xVisium(here("SpatialBenchVisium/preprocessed/FFPE_CytAssist_WT_709/outs/"),
                     type = "sparse",   # use sparse (not HDF5) format
                     data = "raw",     
                     images = "lowres", # specify which image(s) to include
                     load = TRUE)      # specify whether or not to load image(s)


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
    colData(spe)[,"sample_id"] <- rep("FFPE_CytAssist_WT_709",dim(spe)[2])
   
    spe3 <- spe    

    overlap_genes <- intersect(intersect(rownames(spe1),rownames(spe2)),intersect(rownames(spe3),rownames(spe1)))

    spe1 <- spe1[which(rownames(spe1) %in% overlap_genes),]
    spe2 <- spe2[which(rownames(spe2) %in% overlap_genes),]
    spe3 <- spe3[which(rownames(spe3) %in% overlap_genes),]

    spe <- cbind(spe1,spe2,spe3)

    rowData(spe)[,"gene_id"] <- rownames(spe)
    rowData(spe)[,"gene_name"] <- rowData(spe)[,"symbol"]

    spe <- logNormCounts(spe)

    colData(spe)[,"key"] <- paste0(colData(spe)[,"sample_id"], "_", colnames(spe))

    save(spe, file=here("SpatialBenchVisium_probe_based","results", "SpatialBenchVisium_probe_based_spe_qc.Rdata"))
    """
}
