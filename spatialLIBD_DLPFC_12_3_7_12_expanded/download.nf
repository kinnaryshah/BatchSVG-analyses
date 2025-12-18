#!/usr/bin/env nextflow

/*
 * download spatialLIBD DLPFC dataset with 12 samples
 */
process download_data {

    publishDir 'results', mode: 'copy'

    input:
        val data_name

    script:
    """
    #!/usr/bin/env Rscript
    library(SpatialExperiment)
    library(scater)
    library(spatialLIBD)
    library(SpotSweeper)
    library(here)

    spe <- fetch_data("spe")
    colData(spe) <- colData(spe)[,c("sample_id","layer_guess_reordered","sum_umi","sum_gene","subject","position","replicate","subject_position","discard","key","cell_count","in_tissue","array_row","array_col")]
    
    is_mito <- grepl("^mt-", rowData(spe)[,"gene_name"], ignore.case=T)

    # calculate per-spot QC metrics and store in colData
    spe <- addPerCellQC(spe, subsets = list(mito = is_mito))

    # find local outliers with SpotSweeper
    spe <- localOutliers(spe, metric = "sum", direction = "lower", log = TRUE)
    spe <- localOutliers(spe, metric = "detected", direction = "lower", log = TRUE)
    spe <- localOutliers(spe, metric = "subsets_mito_percent", direction = "higher", log = FALSE)

    colData(spe)[,"local_outliers"] = colData(spe)[,"sum_outliers"] | colData(spe)[,"detected_outliers"] | colData(spe)[,"subsets_mito_percent_outliers"]

    spe <- spe[,colData(spe)[,"local_outliers"]==FALSE]

    spe <- spe[,colData(spe)[,"sample_id"] %in% c("151509", "151671", "151676")]

    save(spe, file=here("$data_name","results", "${data_name}_spe_qc.Rdata"))
    """
}