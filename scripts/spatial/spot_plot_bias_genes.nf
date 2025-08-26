#!/usr/bin/env nextflow

/*
 * plot biased genes
 */
process spot_plot_bias {

    publishDir 'results', mode: 'copy'

    input:
        val data_name
        val bias_gene

    script:
    """
    #!/usr/bin/env Rscript
    library(ggspavis)
    library(here)
    library(ggplot2)

    print(paste0("gene: ",$bias_gene))
    
    load(file=here("$data_name","results", "${data_name}_spe_qc.Rdata"))
    print(dim(spe))

    # something is wrong when i use header: true in creating the channel
    if ($bias_gene != "gene_id"){

        splot <- plotCoords(spe, annotate=$bias_gene, assay="logcounts", 
          sample_id="sample_id", point_size=.1)

        pdf(file=here("$data_name","plots",paste0("spot_plot_${data_name}_",$bias_gene,".pdf"))) # issue with quotation marks around bias_gene val 

        plot(splot)

        dev.off()
    }

    """
}
