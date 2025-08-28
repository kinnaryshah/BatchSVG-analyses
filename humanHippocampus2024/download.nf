#!/usr/bin/env nextflow

/*
 * download HPC dataset from ExperimentHub package
 */
process download_data {

    publishDir 'results', mode: 'copy'

    input:
        val data_name

    script:
    """
    #!/usr/bin/env Rscript
    library(ExperimentHub)
    library(here)
    library(dplyr)

    ehub <- ExperimentHub()
    spe <- ehub[["EH9605"]]

    fix_order <- distinct(
        as.data.frame(colData(spe)), slide, array, brnum, sample_id, 
        position, sex) %>% 
        arrange(slide, array)
    sub4 <- fix_order[,"sample_id"][c(14,16,20,21)]
    spe <- spe[,colData(spe)[,"sample_id"] %in% sub4]
    rownames(spe) <- rowData(spe)[,"gene_id"] # need for ggspavis

    print(dim(spe))

    save(spe, file=here("$data_name","results", "${data_name}_spe.Rdata"))
    """
}
