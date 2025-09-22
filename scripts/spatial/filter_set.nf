#!/usr/bin/env nextflow

/*
 * filter out biased genes from SVG list
 */
process filter_out {

    publishDir 'results', mode: 'copy'

    input:
        val data_name
        val batch
        val dev_thres
        val rank_thres

    script:
    """
    #!/usr/bin/env Rscript
    library(BatchSVG)
    library(here)
    
    tab <- read.csv(file=here("$data_name","results", "${data_name}_${batch}_${dev_thres}_${rank_thres}_bias_genes.csv"))

    # load SVGs
    svgs <- read.csv(here("$data_name","results","${data_name}_svgs.csv"), row.names=1,
        check.names=F)
    print(nrow(svgs))

    svgs_filt <- setdiff(svgs[,"gene_id"], tab[,"gene_id"])
    svgs <- svgs[svgs[,"gene_id"] %in% svgs_filt, ]
    print(nrow(svgs))

    write.csv(svgs, file=here("$data_name","results", "${data_name}_${batch}_${dev_thres}_${rank_thres}_filt_svgs.csv"), row.names=F)

    """
}

process filter_out_csv {

    publishDir 'results', mode: 'copy'

    input:
        val data_name
        val bias_csv_concat
        

    script:
    """
    #!/usr/bin/env Rscript
    library(BatchSVG)
    library(here)
    
    tab <- read.csv(file=here("$data_name","results", "$bias_csv_concat"))

    # load SVGs
    svgs <- read.csv(here("$data_name","results","${data_name}_svgs.csv"), row.names=1,
        check.names=F)
    print(nrow(svgs))

    svgs_filt <- setdiff(svgs[,"gene_id"], tab[,"gene_id"])
    svgs <- svgs[svgs[,"gene_id"] %in% svgs_filt, ]
    print(nrow(svgs))

    write.csv(svgs, file=here("$data_name","results", "${data_name}_concat_filt_svgs.csv"), row.names=F)

    """
}
