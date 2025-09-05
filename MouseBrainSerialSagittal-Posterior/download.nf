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

    spe <- read10xVisium(here("MouseBrainSerialSagittal-Posterior/preprocessed/section1/outs/"),
                     type = "sparse",   # use sparse (not HDF5) format
                     data = "raw",     
                     images = "lowres", # specify which image(s) to include
                     load = TRUE)      # specify whether or not to load image(s)

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

    spe1 <- spe

    spe <- read10xVisium(here("MouseBrainSerialSagittal-Posterior/preprocessed/section2/outs/"),
                     type = "sparse",   # use sparse (not HDF5) format
                     data = "raw",     
                     images = "lowres", # specify which image(s) to include
                     load = TRUE)      # specify whether or not to load image(s)

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

    spe2 <- spe

    overlap_genes <- intersect(rownames(spe1),rownames(spe2))

    spe1 <- spe1[which(rownames(spe1) %in% overlap_genes),]
    spe2 <- spe2[which(rownames(spe2) %in% overlap_genes),]

    colData(spe2)[,"sample_id"] <- rep("sample02",dim(spe2)[2])

    spe <- cbind(spe1,spe2)

    save(spe, file=here("$data_name","results", "${data_name}_spe_qc.Rdata"))
    """
}
