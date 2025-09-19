#!/usr/bin/env nextflow

/*
 * Pipeline parameters
 */
params.data_name = 'SpatialBenchVisium_708_709_713'
params.batch1 = "sample_id"
params.batch2 = "preservation"
params.batch3 = "experiment"
params.dev_thres1 = 4
params.rank_thres1 = 7
params.dev_thres2 = 4
params.rank_thres2 = 7
params.dev_thres3 = 4
params.rank_thres3 = 7
params.bias_csv1 = "./results/SpatialBenchVisium_708_709_713_sample_id_4_7_bias_genes.csv"
params.bias_csv2 = "./results/SpatialBenchVisium_708_709_713_preservation_4_7_bias_genes.csv"
params.bias_csv3 = "./results/SpatialBenchVisium_708_709_713_experiment_4_7_bias_genes.csv"
params.bias_csv_concat = "./results/SpatialBenchVisium_708_709_713_concat_bias_genes.csv"


// Include modules
include { download_data } from './download.nf'
include { qc_pipeline } from './../scripts/spatial/qc.nf'
include { feat_sel } from './../scripts/spatial/feature_select.nf'
include { vis_thres } from './../scripts/spatial/vis_thresholds.nf'
include { bias } from './../scripts/spatial/bias_genes.nf'
include { filter_out } from './../scripts/spatial/filter_set.nf'
include { spot_plot_bias } from './../scripts/spatial/bias_genes_plot.nf'
include { bar_plot_bias } from './../scripts/spatial/bias_genes_plot.nf'


workflow {

    // assume I've run 4 download .sh scripts

    download_data(params.data_name)

    // run nnSVG shell script

    // feat_sel(params.data_name, params.batch1)

    // vis_thres(params.data_name, params.batch1, params.dev_thres1, params.rank_thres1)

    // bias(params.data_name, params.batch1, params.dev_thres1, params.rank_thres1)

    // feat_sel(params.data_name, params.batch2)

    // vis_thres(params.data_name, params.batch2, params.dev_thres2, params.rank_thres2)

    // bias(params.data_name, params.batch2, params.dev_thres2, params.rank_thres2)   

    // feat_sel(params.data_name, params.batch3)

    // vis_thres(params.data_name, params.batch3, params.dev_thres3, params.rank_thres3)

    // bias(params.data_name, params.batch3, params.dev_thres3, params.rank_thres3)

    // concat the lists together to create one BatchSVG biased genes csv file
    // concat(params.data_name, params.bias_csv1, params.bias_csv2, params.bias_csv3)

    // from nextflow tutorial: create a channel for inputs from a CSV file
    // bias_gene_ch = Channel.fromPath(params.bias_csv_concat)
    //                    .splitCsv()
    //                    .map { line -> line[0] }

    // spot_plot_bias(params.data_name, bias_gene_ch)
    // bar_plot_bias(params.data_name, bias_gene_ch)

    // filter_out(params.data_name, params.batch, params.dev_thres, params.rank_thres)

}
