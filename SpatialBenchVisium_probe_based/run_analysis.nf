#!/usr/bin/env nextflow

/*
 * Pipeline parameters
 */
params.data_name = 'SpatialBenchVisium_probe_based'
params.batch = "sample_id"
params.dev_thres = 4
params.rank_thres = 7
params.bias_csv = "./results/SpatialBenchVisium_probe_based_sample_id_4_7_bias_genes.csv"

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

    // download_data(params.data_name)

    // run nnSVG shell script

    // feat_sel(params.data_name, params.batch)

    // vis_thres(params.data_name, params.batch, params.dev_thres, params.rank_thres)

    // bias(params.data_name, params.batch, params.dev_thres, params.rank_thres)

    // from nextflow tutorial: create a channel for inputs from a CSV file
    // bias_gene_ch = Channel.fromPath(params.bias_csv)
    //                    .splitCsv()
    //                    .map { line -> line[0] }

    // spot_plot_bias(params.data_name, bias_gene_ch)
    // bar_plot_bias(params.data_name, bias_gene_ch)

    filter_out(params.data_name, params.batch, params.dev_thres, params.rank_thres)

}
