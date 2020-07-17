Analysis of CARNIVAL results
================
Alberto Valdeolivas: <alberto.valdeolivas@bioquant.uni-heidelberg.de>;
Date:
16/07/2020

### License Info

This program is free software: you can redistribute it and/or modify it
under the terms of the GNU General Public License as published by the
Free Software Foundation, either version 3 of the License, or (at your
option) any later version.

This program is distributed in the hope that it will be useful, but
WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General
Public License for more details.

Please check <http://www.gnu.org/licenses/>.

## Introduction

The present script takes the CARNIVAL results inferred from
transcriptomics data from intestinal organoids treated with Sars-CoV-2
for 60 hours and makes some enrichment analysis and clustering.

The original transcriptomic dataset is coming from the following
publication: [Lamers et
al. 2020](https://science.sciencemag.org/content/369/6499/50).

The differential expression analysis was conducted by Martina Poletti
(<Martina.Poletti@earlham.ac.uk>)

## Getting Ready

### Libraries and functions

We first load the required libraries and some important functions to
analyse CARNiVAL output.

``` r
library(readr)
library(fgsea)
library(dplyr)
library(ggplot2)
library(piano)
library(tidyr)
library(ggsci)

## Function to extract the nodes that appear in CARNIVAL network and the 
## background genes (all genes present in the prior knowledge network).
## It returns a list with two objects: the success and the background genes.
extractCARNIVALnodes <- function(CarnivalResults){

    CarnivalNetwork <- 
        as.data.frame(CarnivalResults$weightedSIF, stringsAsFactors = FALSE)
    
    colnames(CarnivalNetwork) <- c("source", "sign", "target", "Weight")

    ## We define the set of nodes interesting for our condition
    sucesses <- unique(c(gsub("_.*","",CarnivalNetwork$source), 
        gsub("_.*","",CarnivalNetwork$target)))

    CarnivalAttributes <- as.data.frame(CarnivalResults$nodesAttributes, 
        stringsAsFactors = FALSE)

    ## We define the background as all the genes in our prior knowledge network.
    bg <- unique(gsub("_.*","",CarnivalAttributes$Node))     
    
    return(list(sucesses = sucesses, bg= bg))
}

### Function to print a barplot with the enriched pathways.
BarplotEnrichment <- function(PathwaysSelect, Interesting_pathways){ 
    
    p <- ggplot(PathwaysSelect, aes(x = reorder(pathway, pvalue), 
            y = -log10(pvalue))) + 
        geom_bar(aes(fill = mean_stat), stat = "identity") +
        scale_fill_gradient2(low = "darkblue", high = "indianred", 
            mid = "whitesmoke", midpoint = 0) + 
        theme_minimal() +
        theme(axis.text.x = element_text(angle = 90, hjust = 1, 
            colour = ifelse(levels(reorder(PathwaysSelect$pathway, 
                PathwaysSelect$pvalue)) %in% Interesting_pathways, 
                "red", "grey40"),
            face = ifelse(levels(reorder(PathwaysSelect$pathway, 
                PathwaysSelect$pvalue)) %in% Interesting_pathways, 
                "bold", "plain")),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank()) + 
        xlab("")
    return(p)
}    
```

### Reading CARNIVAL results

We then read the CARNIVAL results generated in the previous script. We
define two different gene sets in order tor conduct the enrichment. The
first set contains the nodes that appear in the CARNIVAL output and are
therefore relevant in the context of our input transcriptomic data. The
second set contains all the genes in our prior knowledge network which
are used as the backgroud.

``` r
carnival_results_inhibition <- 
    readRDS("Carnival_Results/carnival_results_top50tf_pleitropic_minsize15.rds")
nodes_carnival_inhibition <- extractCARNIVALnodes(carnival_results_inhibition)

carnival_results_activation <- 
    readRDS("Carnival_Results/carnival_results_top50tf_pleitropic_minsize15_activation.rds")
nodes_carnival_activation <- extractCARNIVALnodes(carnival_results_activation)
```

### Reading Pathway data sets from MSigDB

We downloaded from MSigDB <https://www.gsea-msigdb.org/> the following
dataset: c2.cp.v7.0.symbols.gmt. It contains several pathways from
different resources and the genes that are known to be involved in those
pathways.

``` r
Pathway_signatures <- loadGSC("InputData/c2.cp.v7.0.symbols.gmt")
```

### Reading and formatting statistic from DEG

We read the results from the differential expression analysis. The
statistic of the genes will be mapped later on in the different
significant pathways. Maybe, this is not very informative in this
context. I am going to also check with the Average node activity from
CARNIVAL.

``` r
DEA_results <- read_csv("InputData/Clevers_degs_all.csv") %>%
   dplyr::select(Gene, stat)
```

    ## Parsed with column specification:
    ## cols(
    ##   Gene = col_character(),
    ##   baseMean = col_double(),
    ##   log2FoldChange = col_double(),
    ##   lfcSE = col_double(),
    ##   stat = col_double(),
    ##   pvalue = col_double(),
    ##   padj = col_double()
    ## )

``` r
mean_stat <- unlist(lapply(Pathway_signatures$gsc, function(x, DEA_results) {
    genes_matching <- x[which(x %in% DEA_results$Gene)]
    mean_genes <- dplyr::filter(DEA_results, Gene %in% genes_matching) %>%
        dplyr::pull(stat) %>% mean(na.rm = TRUE)
    return(mean_genes)
}, DEA_results = DEA_results)) %>% 
    as.data.frame() %>% 
    tibble::rownames_to_column(var = "pathway") %>%
    dplyr::rename(mean_stat = ".")

### Trying with Carnival node activity
carnival_avg_node_act <- carnival_results_inhibition$nodesAttributes %>% 
    as.data.frame() %>%
    dplyr::mutate(AvgAct = as.numeric(AvgAct))

mean_stat_carni <- 
    unlist(lapply(Pathway_signatures$gsc, function(x, carnival_avg_node_act) {
    genes_matching <- x[which(x %in% carnival_avg_node_act$Node)]
    mean_genes <- dplyr::filter(carnival_avg_node_act, Node %in% genes_matching) %>%
        dplyr::pull(AvgAct) %>% mean(na.rm = TRUE)
    return(mean_genes)
}, carnival_avg_node_act = carnival_avg_node_act)) %>% 
    as.data.frame() %>% 
    tibble::rownames_to_column(var = "pathway") %>%
    dplyr::rename(mean_stat = ".")
```

### Performing Enrichment Analysis and plotting the Results

Using the Piano R package, we run a gene set analysis (GSA) based on a
list of significant genes (CARNIVAL nodes) and a gene set collection
(background). It uses Fisher’s exact test.

First for CARNIVAL results when we assume that the viral-host proteins
interactions are inhibitory:

``` r
hyper_results_inhibition <- 
    runGSAhyper(genes = nodes_carnival_inhibition$sucesses, 
        # pvalues = rep(0,length(nodes_carnival_inhibition$sucesses)),
        pcutoff = 0.01, 
        universe = nodes_carnival_inhibition$bg, gsc = Pathway_signatures,
        gsSizeLim = c(4,Inf), adjMethod = "fdr")
```

    ## Warning in runGSAhyper(genes = nodes_carnival_inhibition$sucesses, pcutoff
    ## = 0.01, : there are genes in gsc that are not in the universe, these will be
    ## removed before analysis

    ## Analyzing the overrepresentation of 85 genes of interest in 2054 gene sets, using a background of 6811 non-interesting genes.

``` r
enriched_pathways_inhibition <- hyper_results_inhibition$resTab %>% 
    as.data.frame() %>% 
    tibble::rownames_to_column(var = "pathway") %>% 
    as_tibble() %>% 
    dplyr::filter(`Adjusted p-value` <= 0.01) %>% 
    dplyr::inner_join(mean_stat_carni) %>% 
    dplyr::select(pathway, `p-value`, `Adjusted p-value`, mean_stat) %>%
    dplyr::rename(pvalue = `p-value`, AdjPvalu = `Adjusted p-value`) %>% 
    dplyr::mutate(pathway = as.factor(pathway))
```

    ## Joining, by = "pathway"

``` r
interesting_pathways_inhibition <- c()

p_inhibition <- BarplotEnrichment(enriched_pathways_inhibition, 
    interesting_pathways_inhibition)
```

![](02_analysis_carnival_results_files/figure-gfm/unnamed-chunk-7-1.png)<!-- -->

Then, for CARNIVAL when we assume that the viral-host proteins
interactions are stimulations:

``` r
hyper_results_activation <- 
    runGSAhyper(genes = nodes_carnival_activation$sucesses, 
        # pvalues = rep(0,length(nodes_carnival_activation$sucesses)),
        pcutoff = 0.01, 
        universe = nodes_carnival_activation$bg, gsc = Pathway_signatures,
        gsSizeLim = c(4,Inf), adjMethod = "fdr")
```

    ## Warning in runGSAhyper(genes = nodes_carnival_activation$sucesses, pcutoff
    ## = 0.01, : there are genes in gsc that are not in the universe, these will be
    ## removed before analysis

    ## Analyzing the overrepresentation of 86 genes of interest in 2054 gene sets, using a background of 6810 non-interesting genes.

``` r
enriched_pathways_activation <- hyper_results_activation$resTab %>% 
    as.data.frame() %>% 
    tibble::rownames_to_column(var = "pathway") %>% 
    as_tibble() %>% 
    dplyr::filter(`Adjusted p-value` <= 0.01) %>% 
    dplyr::inner_join(mean_stat_carni) %>% 
    dplyr::select(pathway, `p-value`, `Adjusted p-value`, mean_stat) %>%
    dplyr::rename(pvalue = `p-value`, AdjPvalu = `Adjusted p-value`) %>% 
    dplyr::mutate(pathway = as.factor(pathway))
```

    ## Joining, by = "pathway"

``` r
interesting_pathways_activation <- c()

p_activation <- BarplotEnrichment(enriched_pathways_activation, 
    interesting_pathways_activation)
```

![](02_analysis_carnival_results_files/figure-gfm/unnamed-chunk-10-1.png)<!-- -->

### Comparing Enrichment between activation and inhibition

Now I compare the enriched pathways that change the most between both
assumptions: host-viral interactions are inhibitory or stimulatory.

``` r
results_activation <-   hyper_results_activation$resTab %>% 
    as.data.frame() %>% 
    tibble::rownames_to_column(var = "pathway") %>% 
    as_tibble() %>% 
    dplyr::mutate(LogpValue_activation = -log(`p-value`), 
        p_value_activation = `p-value`) %>% 
    dplyr::select(pathway, p_value_activation, LogpValue_activation) 

results_inhibition <-   hyper_results_inhibition$resTab %>% 
    as.data.frame() %>% 
    tibble::rownames_to_column(var = "pathway") %>% 
    as_tibble() %>% 
    dplyr::mutate(LogpValue_inhibition = -log(`p-value`), 
        p_value_inhibition = `p-value`) %>% 
    dplyr::select(pathway, p_value_inhibition, LogpValue_inhibition) 

results_join <- 
    dplyr::inner_join(results_activation, results_inhibition) %>% 
    dplyr::filter(p_value_activation < 0.05 | p_value_inhibition < 0.05) %>% 
    dplyr::mutate(diff = LogpValue_activation - LogpValue_inhibition) %>% 
    dplyr::top_n(15, wt = abs(diff)) %>% 
    dplyr::arrange(desc(abs(diff)))
```

    ## Joining, by = "pathway"

``` r
results_join_longer <- results_join %>% 
    dplyr::select(pathway, diff, LogpValue_activation, LogpValue_inhibition) %>% 
    dplyr::rename(activation = "LogpValue_activation", 
        inhibition = "LogpValue_inhibition") %>% 
    pivot_longer(-c(pathway,diff), values_to = "LogPvalue", 
        names_to = "Interaction")  
    
point_plot_activation_vs_inhibtion <- 
    ggplot(results_join_longer, aes(reorder(pathway, abs(diff)), LogPvalue)) + 
    geom_point(aes(color = Interaction), size = 3) + 
    coord_flip() + 
    theme_minimal() + 
    theme(legend.position = "bottom",  legend.justification = "center") +
    scale_color_lancet() +
    theme(axis.text.x = element_text(hjust = 1, size =8.5),
        axis.text.y = element_text(size = 7),
        panel.grid.minor = element_blank()) + 
    xlab("Pathways") + ylab("-Log (p-value)") +
    geom_hline(yintercept = -log(0.05), linetype="dashed", 
                color = "#2F4F4F", size=0.5)
```

``` r
point_plot_activation_vs_inhibtion
```

![](02_analysis_carnival_results_files/figure-gfm/unnamed-chunk-12-1.png)<!-- -->

## Session Info Details

    ## R version 4.0.1 (2020-06-06)
    ## Platform: x86_64-pc-linux-gnu (64-bit)
    ## Running under: Ubuntu 19.10
    ## 
    ## Matrix products: default
    ## BLAS:   /usr/lib/x86_64-linux-gnu/openblas/libblas.so.3
    ## LAPACK: /usr/lib/x86_64-linux-gnu/libopenblasp-r0.3.7.so
    ## 
    ## locale:
    ##  [1] LC_CTYPE=en_GB.UTF-8       LC_NUMERIC=C              
    ##  [3] LC_TIME=en_GB.UTF-8        LC_COLLATE=en_GB.UTF-8    
    ##  [5] LC_MONETARY=en_GB.UTF-8    LC_MESSAGES=en_GB.UTF-8   
    ##  [7] LC_PAPER=en_GB.UTF-8       LC_NAME=C                 
    ##  [9] LC_ADDRESS=C               LC_TELEPHONE=C            
    ## [11] LC_MEASUREMENT=en_GB.UTF-8 LC_IDENTIFICATION=C       
    ## 
    ## attached base packages:
    ## [1] stats     graphics  grDevices utils     datasets  methods   base     
    ## 
    ## other attached packages:
    ## [1] ggsci_2.9     tidyr_1.1.0   piano_2.4.0   ggplot2_3.3.1 dplyr_1.0.0  
    ## [6] fgsea_1.14.0  readr_1.3.1  
    ## 
    ## loaded via a namespace (and not attached):
    ##  [1] Rcpp_1.0.4.6         lattice_0.20-41      visNetwork_2.0.9    
    ##  [4] relations_0.6-9      gtools_3.8.2         digest_0.6.25       
    ##  [7] mime_0.9             slam_0.1-47          R6_2.4.1            
    ## [10] evaluate_0.14        pillar_1.4.4         gplots_3.0.3        
    ## [13] rlang_0.4.6          data.table_1.12.8    gdata_2.18.0        
    ## [16] DT_0.13              Matrix_1.2-18        rmarkdown_2.2       
    ## [19] labeling_0.3         shinyjs_1.1          sets_1.0-18         
    ## [22] BiocParallel_1.22.0  stringr_1.4.0        htmlwidgets_1.5.1   
    ## [25] igraph_1.2.5         munsell_0.5.0        shiny_1.4.0.2       
    ## [28] compiler_4.0.1       httpuv_1.5.4         xfun_0.14           
    ## [31] pkgconfig_2.0.3      BiocGenerics_0.34.0  marray_1.66.0       
    ## [34] htmltools_0.4.0      tidyselect_1.1.0     tibble_3.0.1        
    ## [37] gridExtra_2.3        crayon_1.3.4         withr_2.2.0         
    ## [40] later_1.1.0.1        bitops_1.0-6         grid_4.0.1          
    ## [43] jsonlite_1.6.1       xtable_1.8-4         gtable_0.3.0        
    ## [46] lifecycle_0.2.0      magrittr_1.5         scales_1.1.1        
    ## [49] KernSmooth_2.23-17   stringi_1.4.6        farver_2.0.3        
    ## [52] promises_1.1.1       limma_3.44.3         ellipsis_0.3.1      
    ## [55] generics_0.0.2       vctrs_0.3.1          fastmatch_1.1-0     
    ## [58] tools_4.0.1          Biobase_2.48.0       glue_1.4.1          
    ## [61] purrr_0.3.4          hms_0.5.3            parallel_4.0.1      
    ## [64] fastmap_1.0.1        yaml_2.2.1           colorspace_1.4-1    
    ## [67] cluster_2.1.0        caTools_1.18.0       shinydashboard_0.7.1
    ## [70] knitr_1.28
