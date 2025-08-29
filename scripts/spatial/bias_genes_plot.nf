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

process bar_plot_bias {

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
    library(SpatialExperiment)
    library(tidyverse)

    print(paste0("gene: ",$bias_gene))
    
    load(file=here("$data_name","results", "${data_name}_spe_qc.Rdata"))
    print(dim(spe))

    # something is wrong when i use header: true in creating the channel
    if ($bias_gene != "gene_id"){

        nonzero_list <- c()
        expr_list <- c()
        for(i in unique(colData(spe)[,"sample_id"])){
    
            spe_sub <- spe[rownames(spe)==$bias_gene,colData(spe)[,"sample_id"] == i]
            nonzero <- sum(logcounts(spe_sub) > 0) / dim(spe_sub)[2]
            nonzero_list <- append(nonzero_list, nonzero)
    
            expr <- mean(logcounts(spe_sub))
            expr_list <- append(expr_list, expr)
    
    }
        
        df <- data.frame(sample = unique(colData(spe)[,"sample_id"]),
                        nonzero_frac = nonzero_list,
                        avg_expr = expr_list)

        print(head(df))

        # https://stackoverflow.com/questions/64298884/creating-barplot-with-dual-y-axis
        #Scaling factor
        sf <- max(df[,"nonzero_frac"])/max(df[,"avg_expr"])
        #Transform
        DF_long <- df %>%
        mutate(avg_expr = avg_expr*sf) %>%
        pivot_longer(names_to = "y_new", values_to = "val", nonzero_frac:avg_expr)

        plot <- ggplot(DF_long, aes(x=sample)) +
        geom_bar( aes(y = val, fill = y_new, group = y_new),
                    stat="identity", position=position_dodge(),
                    color="black", alpha=.6)  +
        scale_fill_manual(values = c("blue", "red")) +
        scale_y_continuous(name = "nonzero_frac",labels = scales::comma,sec.axis = sec_axis(~./sf, name="avg_expr",
                                                                                    labels = scales::comma))+
        labs(fill='variable')+
        theme_bw()+
        theme(legend.position = 'top',
                plot.title = element_text(color='black',face='bold',hjust=0.5),
                axis.text = element_text(color='black',face='bold'),
                axis.title = element_text(color='black',face='bold'),
                legend.text = element_text(color='black',face='bold'),
                legend.title = element_text(color='black',face='bold')) +
            ggtitle($bias_gene) +
            theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

        pdf(file=here("$data_name","plots",paste0("bar_plot_${data_name}_",$bias_gene,".pdf"))) # issue with quotation marks around bias_gene val 

        plot(plot)

        dev.off()

    }

    """
}
