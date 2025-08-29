#!/usr/bin/env nextflow

/*
 * Pipeline parameters
 */
params.data_name = 'spatialLIBD_DLPFC_12'
params.batch = "subject"
params.dev_thres = 4
params.rank_thres = 7
params.bias_csv = "./results/spatialLIBD_DLPFC_12_subject_4_7_bias_genes.csv"

// Include modules
include { download_data } from './download.nf'
include { qc_pipeline } from './../scripts/spatial/qc.nf'
include { feat_sel } from './../scripts/spatial/feature_select.nf'
include { vis_thres } from './../scripts/spatial/vis_thresholds.nf'
include { bias } from './../scripts/spatial/bias_genes.nf'
include { filter_out } from './../scripts/spatial/filter_set.nf'
include { spot_plot_bias } from './../scripts/spatial/bias_genes_plot.nf'


workflow {

    download_data(params.data_name)

    qc_pipeline(params.data_name)

    // normally would have step to run nnSVG, but already have results

    feat_sel(params.data_name, params.batch)

    vis_thres(params.data_name, params.batch, params.dev_thres, params.rank_thres)

    bias(params.data_name, params.batch, params.dev_thres, params.rank_thres)

    // from nextflow tutorial: create a channel for inputs from a CSV file
    bias_gene_ch = Channel.fromPath(params.bias_csv)
                        .splitCsv()
                        .map { line -> line[0] }

    spot_plot_bias(params.data_name, bias_gene_ch)

    filter_out(params.data_name, params.batch, params.dev_thres, params.rank_thres)

}
