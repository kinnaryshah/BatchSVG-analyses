#!/usr/bin/env nextflow

/*
 * run Batch_SVG::feature_select()
 */
process feat_sel {

    publishDir 'results', mode: 'copy'

    input:
        val data_name
        val batch

    script:
    """
    #!/usr/bin/env Rscript
    library(BatchSVG)
    library(here)

    load(file=here("$data_name","results", "${data_name}_spe_qc.Rdata"))
    print(dim(spe))

    # load SVGs
    svgs <- read.csv(here("$data_name","results","${data_name}_svgs.csv"), row.names=1,
        check.names=F)
    print(nrow(svgs))

    # run function
    list_batch_df <- featureSelect(input = spe, 
        batch_effect = "$batch", VGs = svgs[,"gene_id"])
    head(list_batch_df[["$batch"]])

    save(list_batch_df, file=here("$data_name","results", "${data_name}_${batch}_feat_sel.RData"))
    """
}
