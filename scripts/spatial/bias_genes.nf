#!/usr/bin/env nextflow

/*
 * find biased genes
 */
process bias {

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
    
    bias_both <- biasDetect(list_batch_df = list_batch_df, threshold = "both",
        nSD_dev = $dev_thres, nSD_rank = $rank_thres, plot_point_shape = 23, plot_palette = "RdPu",
        plot_text_size = 4)

    tab <- bias_both[["$batch"]][["Table"]][,c("gene_id","gene_name","nSD_bin_dev", "dev_outlier", "nSD_bin_rank", "rank_outlier")]
    write.csv(tab, file=here("$data_name","results", "${data_name}_${batch}_${dev_thres}_${rank_thres}_bias_genes.csv"), row.names=F)

    pdf(here("$data_name","plots", "bias_features_${data_name}_${batch}_${dev_thres}_${rank_thres}.pdf"))
    
    bias_both[["$batch"]][["Plot"]]

    dev.off()

    """
}
