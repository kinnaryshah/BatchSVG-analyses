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

    spe <- read10xVisium(here("SpatialBenchVisium_OCT_manual/preprocessed/OCT_WT_544/outs/"),
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
    colData(spe)[,"sample_id"] <- rep("OCT_WT_544",dim(spe)[2])
    colData(spe)[,"preservation"] <- rep("OCT",dim(spe)[2])
    colData(spe)[,"placement"] <- rep("manual",dim(spe)[2])
    colData(spe)[,"experiment"] <- rep("polyA",dim(spe)[2])
       
    spe1 <- spe

    spe <- read10xVisium(here("SpatialBenchVisium_OCT_manual/preprocessed/OCT_WT_545/outs/"),
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
    colData(spe)[,"sample_id"] <- rep("OCT_WT_545",dim(spe)[2])
    colData(spe)[,"preservation"] <- rep("OCT",dim(spe)[2])
    colData(spe)[,"placement"] <- rep("manual",dim(spe)[2])
    colData(spe)[,"experiment"] <- rep("polyA",dim(spe)[2])
       
    spe2 <- spe    

    spe <- read10xVisium(here("SpatialBenchVisium_OCT_manual/preprocessed/OCT_WT_708/outs/"),
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
    colData(spe)[,"sample_id"] <- rep("OCT_WT_708",dim(spe)[2])
    colData(spe)[,"preservation"] <- rep("OCT",dim(spe)[2])
    colData(spe)[,"placement"] <- rep("manual",dim(spe)[2])
     colData(spe)[,"experiment"] <- rep("polA",dim(spe)[2])
      
    spe3 <- spe   

    spe <- read10xVisium(here("SpatialBenchVisium_OCT_manual/preprocessed/OCT_WT_709/outs/"),
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
    colData(spe)[,"sample_id"] <- rep("OCT_WT_709",dim(spe)[2])
    colData(spe)[,"preservation"] <- rep("OCT",dim(spe)[2])
    colData(spe)[,"placement"] <- rep("manual",dim(spe)[2])
     colData(spe)[,"experiment"] <- rep("polyA",dim(spe)[2])
      
    spe4 <- spe   

    overlap_genes <- Reduce(intersect,  list(v1 = rownames(spe1), 
                                            v2 = rownames(spe2), 
                                            v3 = rownames(spe3),
                                            v4 = rownames(spe4))) 

    print(length(overlap_genes))

    spe1 <- spe1[which(rownames(spe1) %in% overlap_genes),]
    spe2 <- spe2[which(rownames(spe2) %in% overlap_genes),]
    spe3 <- spe3[which(rownames(spe3) %in% overlap_genes),]
    spe4 <- spe4[which(rownames(spe4) %in% overlap_genes),]

    spe <- cbind(spe1,spe2,spe3,spe4)

    rowData(spe)[,"gene_id"] <- rownames(spe)
    rowData(spe)[,"gene_name"] <- rowData(spe)[,"symbol"]

    spe <- logNormCounts(spe)

    colData(spe)[,"key"] <- paste0(colData(spe)[,"sample_id"], "_", colnames(spe))

    save(spe, file=here("SpatialBenchVisium_OCT_manual","results", "SpatialBenchVisium_OCT_manual_spe_qc.Rdata"))
    """
}
