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

    df1 <- read.csv(bias_csv1)
    df2 <- read.csv(bias_csv2)
    df3 <- read.csv(bias_csv3)

    df <- rbind(df1,df2,df3)

    print(dim(df))

    write.csv(df, file=here("$data_name","results", "${data_name}_concat_bias_genes.csv"), row.names=F)
    

    """
}
