#!/usr/bin/env nextflow

/*
 * run quality control on dataset
 */
process qc_pipeline {

    publishDir 'results', mode: 'copy'

    input:
        val data_name

    script:
    """
    #!/usr/bin/env Rscript
    library(spatialLIBD)
    library(here)

    load(file=here("$data_name","results", "${data_name}_spe.Rdata"))

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

    save(spe, file=here("$data_name","results", "${data_name}_spe_qc.Rdata"))
    """
}
