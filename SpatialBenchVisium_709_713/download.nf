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

    spe <- read10xVisium(here("SpatialBenchVisium_709_713/preprocessed/OCT_CytAssist_WT_709/outs/"),
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
    colData(spe)[,"preservation"] <- rep("OCT",dim(spe)[2])
    colData(spe)[,"placement"] <- rep("CytAssist",dim(spe)[2])
    colData(spe)[,"experiment"] <- rep("probe",dim(spe)[2])
       
    spe1 <- spe

    spe <- read10xVisium(here("SpatialBenchVisium_709_713/preprocessed/FFPE_2_WT_709/outs/"),
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
    colData(spe)[,"preservation"] <- rep("FFPE",dim(spe)[2])
    colData(spe)[,"placement"] <- rep("manual",dim(spe)[2])
    colData(spe)[,"experiment"] <- rep("probe",dim(spe)[2])
       
    spe2 <- spe

    spe <- read10xVisium(here("SpatialBenchVisium_709_713/preprocessed/FFPE_CytAssist_WT_709/outs/"),
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
    colData(spe)[,"preservation"] <- rep("FFPE",dim(spe)[2])
    colData(spe)[,"placement"] <- rep("CytAssist",dim(spe)[2])
    colData(spe)[,"experiment"] <- rep("probe",dim(spe)[2])
       
    spe3 <- spe    

    spe <- read10xVisium(here("SpatialBenchVisium_709_713/preprocessed/FFPE_2_WT_713/outs/"),
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
    colData(spe)[,"sample_id"] <- rep("FFPE_2_WT_713",dim(spe)[2])
    colData(spe)[,"preservation"] <- rep("FFPE",dim(spe)[2])
    colData(spe)[,"placement"] <- rep("manual",dim(spe)[2])
     colData(spe)[,"experiment"] <- rep("probe",dim(spe)[2])
      
    spe4 <- spe   


    spe <- read10xVisium(here("SpatialBenchVisium_709_713/preprocessed/OCT_CytAssist_WT_713/outs/"),
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
    colData(spe)[,"sample_id"] <- rep("OCT_CytAssist_WT_713",dim(spe)[2])
    colData(spe)[,"preservation"] <- rep("OCT",dim(spe)[2])
    colData(spe)[,"placement"] <- rep("CytAssist",dim(spe)[2])
    colData(spe)[,"experiment"] <- rep("probe",dim(spe)[2])
       
    spe5 <- spe 


    spe <- read10xVisium(here("SpatialBenchVisium_709_713/preprocessed/FFPE_CytAssist_WT_713/outs/"),
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
    colData(spe)[,"sample_id"] <- rep("FFPE_CytAssist_WT_713",dim(spe)[2])
    colData(spe)[,"preservation"] <- rep("FFPE",dim(spe)[2])
    colData(spe)[,"placement"] <- rep("CytAssist",dim(spe)[2])
    colData(spe)[,"experiment"] <- rep("probe",dim(spe)[2])
       
    spe6 <- spe 

    overlap_genes <- Reduce(intersect,  list(v1 = rownames(spe1), 
                                            v2 = rownames(spe2), 
                                            v3 = rownames(spe3),
                                            v4 = rownames(spe4),
                                            v5= rownames(spe5),
                                            v6 = rownames(spe6))) 

    print(length(overlap_genes))

    spe1 <- spe1[which(rownames(spe1) %in% overlap_genes),]
    spe2 <- spe2[which(rownames(spe2) %in% overlap_genes),]
    spe3 <- spe3[which(rownames(spe3) %in% overlap_genes),]
    spe4 <- spe4[which(rownames(spe4) %in% overlap_genes),]
    spe5 <- spe5[which(rownames(spe5) %in% overlap_genes),]
    spe6 <- spe6[which(rownames(spe6) %in% overlap_genes),]

    spe <- cbind(spe1,spe2,spe3,spe4,spe5,spe6)

    rowData(spe)[,"gene_id"] <- rownames(spe)
    rowData(spe)[,"gene_name"] <- rowData(spe)[,"symbol"]

    spe <- logNormCounts(spe)

    colData(spe)[,"key"] <- paste0(colData(spe)[,"sample_id"], "_", colnames(spe))

    save(spe, file=here("SpatialBenchVisium_709_713","results", "SpatialBenchVisium_709_713_spe_qc.Rdata"))
    """
}
