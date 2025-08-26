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

    spe <- fetch_data(type = "spe")

    print(dim(spe))

    save(spe, file=here("$data_name","results", "${data_name}_spe.Rdata"))
    """
}
