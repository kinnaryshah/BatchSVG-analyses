#!/usr/bin/env nextflow

/*
 * concat multiple BatchSVG files together
 */
process concat {

    publishDir 'results', mode: 'copy'

    input:
        val data_name
        val bias_csv1
        val bias_csv2
        val bias_csv3

    script:
    """
    #!/usr/bin/env Rscript
    library(BatchSVG)
    library(here)

    df1 <- read.csv(here("$data_name","results", "${data_name}_sample_id_9_5_bias_genes.csv"))
    df2 <- read.csv(here("$data_name","results", "${data_name}_preservation_7_5_bias_genes.csv"))
    df3 <- read.csv(here("$data_name","results", "${data_name}_experiment_12_5_bias_genes.csv"))

    df <- rbind(df1,df2,df3)

    print(dim(df))

    write.csv(df, file=here("$data_name","results", "${data_name}_concat_bias_genes.csv"), row.names=F)
    

    """
}
