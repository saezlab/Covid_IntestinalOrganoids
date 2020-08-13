# Import libraries
library(tibble)
library(dplyr)
library(readr)
library(progeny)
library(dorothea)
library(CARNIVAL)
library(ggplot2)
library(OmnipathR)
library(biomaRt)
library(grid)
library(gridExtra)

# Create prior knowledge network from omnipath
ia_omnipath <- import_omnipath_interactions() %>% as_tibble()

ia_pwextra <- import_pathwayextra_interactions() %>% as_tibble()

ia_kinaseextra <- import_kinaseextra_interactions() %>% as_tibble()

## We bind the datasets
interactions <- as_tibble(
  bind_rows(
    ia_omnipath %>% mutate(type = 'ppi'),
    ia_pwextra %>% mutate(type = 'ppi'),
    ia_kinaseextra %>% mutate(type = 'ppi')))

signed_directed_interactions <- 
  dplyr::filter(interactions, consensus_direction==1) %>%
  filter(consensus_stimulation == 1 | consensus_inhibition == 1) %>% 
  dplyr::mutate(sign = if_else(consensus_stimulation==1,1,-1))  %>%
  dplyr::select(source_genesymbol, sign,  target_genesymbol) %>%
  dplyr::rename(source ="source_genesymbol", target ="target_genesymbol")

carnival_pkn <- signed_directed_interactions %>%
  dplyr::distinct(source, target, .keep_all = TRUE)

all_source_nodes <- unique(carnival_pkn$source)
all_target_nodes <- unique(carnival_pkn$target)
all_nodes_network <- unique(c(all_source_nodes,all_target_nodes))

# Read in viral mirna - host protein network
host_viral_interactions <- 
  read_tsv(file = "mirna_to_human_unfiltered.txt")

## Filter this networks for human proteins in the prior knowledge network
## We assume that they are inhibitions (sign:-1)
host_viral_interactions_filter <- host_viral_interactions %>% 
  dplyr::filter(human_target_gene %in% all_source_nodes) %>% 
  dplyr::select(viral_miRNA_acc, human_target_gene) %>% 
  dplyr::distinct() %>% 
  dplyr::mutate(sign = -1) %>% 
  dplyr::rename(source = "viral_miRNA_acc", target = "human_target_gene") %>%
  dplyr::select(source, sign, target)

## Prior Knowledge Network
carnival_pkn_hostvirus <- 
  bind_rows(carnival_pkn, host_viral_interactions_filter)

## We load Dorothea Regulons
data(dorothea_hs, package = "dorothea")
regulons <- dorothea_hs %>%
  dplyr::filter(confidence %in% c("A", "B","C"))

## Perturbation (Viral proteins are active)
viral_proteins <- unique(host_viral_interactions_filter$source)
viral_proteins_perturbation <- 
  data.frame(viral_proteins = viral_proteins, sign = 1) %>% 
  tibble::column_to_rownames(var = "viral_proteins") %>%
  t() %>% as.data.frame()

## Top TFs
# read in degs
all_genes_dea <- read_csv("Clevers_degs_all.csv") %>% 
  dplyr::filter(!(Gene %in% c("Mar-02","RPL21")))

all_genes_dea_stat <-  all_genes_dea %>% 
  dplyr::select(Gene, stat) %>% 
  dplyr::filter(!is.na(stat)) %>% 
  column_to_rownames(var = "Gene") %>%
  as.matrix()

# TFs have at least 5 targets
tf_activities_stat_pleitropic <- 
  dorothea::run_viper(all_genes_dea_stat, regulons,
                      options =  list(minsize = 5, eset.filter = FALSE, 
                                  verbose = FALSE, nes = TRUE, pleiotropy= TRUE))

tf_top50_pleitropic <- tf_activities_stat_pleitropic %>%
  as.data.frame() %>% 
  rownames_to_column(var = "TF") %>%
  dplyr::filter(TF %in% all_nodes_network) %>% 
  dplyr::rename(NES = "stat") %>%
  dplyr::top_n(50, wt = abs(NES)) %>%
  dplyr::arrange(NES) %>%
  tibble::column_to_rownames(var = "TF")  %>%
  t() %>% as.data.frame()

## Progeny Weigths
# Changed from top 100 to top 100
pathways_inputCarnival <- 
  t(progeny(as.matrix(all_genes_dea_stat), 
            scale=TRUE, organism="Human", top = 100, perm = 10000, z_scores = FALSE))
colnames(pathways_inputCarnival) <- "Activity"

progeny_weigths <- pathways_inputCarnival %>% t()

# Running Carnival
carnival_results_top50tf_pleitropic_minsize5 <-
  runCARNIVAL(
    measObj=tf_top50_pleitropic,
    inputObj = viral_proteins_perturbation,
    netObj=carnival_pkn_hostvirus,
    weightObj=progeny_weigths,
    timelimit = 7200,
    solverPath="/home/alvaldeolias/Downloads/cplex",
    solver = "cplex"
)

saveRDS(carnival_results_top50tf_pleitropic_minsize5, 
    file = "carnival_results_top50tf_pleitropic_minsize5.rds")

