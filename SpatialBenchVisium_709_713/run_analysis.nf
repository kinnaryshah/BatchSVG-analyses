#!/usr/bin/env nextflow

/*
 * Pipeline parameters
 */
params.data_name = 'SpatialBenchVisium_709_713'
params.batch1 = "sample_id"
params.batch2 = "preservation"
params.batch3 = "placement"
params.dev_thres1 = 9
params.rank_thres1 = 5
params.dev_thres2 = 7
params.rank_thres2 = 5
params.dev_thres3 = 12
params.rank_thres3 = 5
params.bias_csv1 = "./results/SpatialBenchVisium_709_713_sample_id_9_5_bias_genes.csv"
params.bias_csv2 = "./results/SpatialBenchVisium_709_713_preservation_7_5_bias_genes.csv"
params.bias_csv3 = "./results/SpatialBenchVisium_709_713_placement_12_5_bias_genes.csv"
params.bias_csv_concat = "SpatialBenchVisium_709_713_concat_bias_genes.csv"


// Include modules
include { download_data } from './download.nf'
include { qc_pipeline } from './../scripts/spatial/qc.nf'
include { feat_sel } from './../scripts/spatial/feature_select.nf'
include { vis_thres } from './../scripts/spatial/vis_thresholds.nf'
include { bias } from './../scripts/spatial/bias_genes.nf'
include { filter_out } from './../scripts/spatial/filter_set.nf'
include { spot_plot_bias } from './../scripts/spatial/bias_genes_plot.nf'
include { bar_plot_bias } from './../scripts/spatial/bias_genes_plot.nf'
include { concat } from './../scripts/spatial/concat_csv.nf'
include { filter_out_csv } from './../scripts/spatial/filter_set.nf'


workflow {

    // assume I've run 4 download .sh scripts

    // download_data(params.data_name)

    // run nnSVG shell script

    // feat_sel(params.data_name, params.batch1)

    // vis_thres(params.data_name, params.batch1, params.dev_thres1, params.rank_thres1)

    // bias(params.data_name, params.batch1, params.dev_thres1, params.rank_thres1)

    // feat_sel(params.data_name, params.batch2)

    // vis_thres(params.data_name, params.batch2, params.dev_thres2, params.rank_thres2)

    // bias(params.data_name, params.batch2, params.dev_thres2, params.rank_thres2)   

    // feat_sel(params.data_name, params.batch3)

    vis_thres(params.data_name, params.batch3, params.dev_thres3, params.rank_thres3)

    // bias(params.data_name, params.batch3, params.dev_thres3, params.rank_thres3)

    // from nextflow tutorial: create a channel for inputs from a CSV file
    // bias_gene_ch1 = Channel.fromPath(params.bias_csv1)
    //                    .splitCsv()
    //                    .map { line -> line[0] }

    // spot_plot_bias(params.data_name, params.batch1, bias_gene_ch1)
    // bar_plot_bias(params.data_name,  params.batch1, bias_gene_ch1)

    // from nextflow tutorial: create a channel for inputs from a CSV file
    // bias_gene_ch2 = Channel.fromPath(params.bias_csv2)
    //                    .splitCsv()
    //                    .map { line -> line[0] }

    // spot_plot_bias(params.data_name, params.batch2, bias_gene_ch2)
    // bar_plot_bias(params.data_name,  params.batch2, bias_gene_ch2)

    // from nextflow tutorial: create a channel for inputs from a CSV file
    // bias_gene_ch3 = Channel.fromPath(params.bias_csv3)
    //                    .splitCsv()
    //                    .map { line -> line[0] }

    // spot_plot_bias(params.data_name, params.batch3, bias_gene_ch3)
    // bar_plot_bias(params.data_name,  params.batch3, bias_gene_ch3)

    // concat the lists together to create one BatchSVG biased genes csv file
    // concat(params.data_name, params.bias_csv1, params.bias_csv2, params.bias_csv3)

    // filter_out_csv(params.data_name, params.bias_csv_concat)

}
