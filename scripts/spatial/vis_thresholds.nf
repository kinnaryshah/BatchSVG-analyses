#!/usr/bin/env nextflow

/*
 * visualize Batch_SVG::feature_select() output
 */
process vis_thres {

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

    load(file=here("$data_name","results", "${data_name}_${batch}_feat_sel.RData"))
    plots <- svg_nSD(list_batch_df = list_batch_df, 
            sd_interval_dev = 4, sd_interval_rank = 7)

    pdf(here("$data_name","plots", "thresholds_${data_name}_${batch}_${dev_thres}_${rank_thres}.pdf"))
    plots[["$batch"]]
    dev.off()
    """
}
