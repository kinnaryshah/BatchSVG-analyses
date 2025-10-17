#!/usr/bin/env nextflow

/*
 * download spatialLIBD DLPFC dataset with 12 samples
 */
process download_data {

    publishDir 'results', mode: 'copy'

    input:
        val data_name

    script:
    """
    #!/usr/bin/env Rscript
    library(spatialLIBD)
    library(here)

    data_load <- "spatialLIBD_DLPFC_12_jacqui"
    load(file = here(data_load, "results", paste0(data_load, "_spe_qc.Rdata")))

    spe <- spe[,colData(spe)[,"sample_id"] %in% c("151509", "151671", "151676")]

    save(spe, file=here("$data_name","results", "${data_name}_spe.Rdata"))
    """
}
