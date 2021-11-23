#' @param args[1] path to cluster file
#' @param args[2] use biogrid genes for background T/F
#' @param args[3] filter out chronic inf seeds T/F

# setup -------------------------------------------------------------------
library(topGO)
library(org.Hs.eg.db)
library(tidyverse)
source("chronic_inflammation_functions.R")

args <- commandArgs(TRUE)
cluster_df = read.csv(args[1], row.names = 1)

biogrid_index = read.delim("../data/biogrid/biogrid_gene_index.txt", header = F)
biogrid_genes = biogrid_index$V2

meta = cluster_df %>%
  dplyr::select(Disease, 
         Graph, 
         Method,
         Resolution,
         Iterations,
         expandFDR,
         expandCutoff) %>%
  distinct()

disease = meta$Disease
iterations =  meta$Iterations
graph = meta$Graph
FDR = meta$expandFDR
cutoff = meta$expandCutoff


ci_seeds = read.delim("../data/hotnet2_heat_files/chronic_inflammation_gene_shot_pubs_greater10.txt",
                      header = F)

ci_seeds = ci_seeds %>%
  filter(V2 == 1) %>%
  pull(V1)

if(args[3] == TRUE){
  cluster_df = cluster_df %>%
    filter(!Gene %in% ci_seeds)
}

keep_clusters = cluster_df %>%
  group_by(Cluster) %>%
  tally() %>%
  filter(n >=5) %>%
  pull(Cluster)

clusters = unique(cluster_df$Cluster)
clusters = clusters[clusters %in% keep_clusters]

all_results = {}

  for(clust in clusters){
    
    print(clust)
    clust_genes = cluster_df %>%
      filter(Cluster == clust) %>%
      pull(Gene)
    
    if(args[2] == TRUE){
      background = biogrid_genes
    } else {
      background = cluster_df$Gene}
  
    results = findEnrichedGOBP(clust_genes, 
                               background,
                               "org.Hs.eg.db",
                               "entrez",
                               min_size = 5,
                               max_size = 100)
    
    
    results$Cluster = clust
    all_results = rbind(all_results, results)
  }

meta_df = do.call(rbind, 
                  replicate(nrow(all_results), 
                            meta, 
                            simplify = FALSE))

all_results = cbind(all_results, meta_df)

if(args[3] == FALSE){
  outdir = "/mnt/research/compbio/krishnanlab/projects/chronic_inflammation/results/biogrid_noRWR/cluster_GOBP_enrichment/expanded_disease_genes_FDR=FALSE_cutoff=0.01"
}else{
  outdir = "/mnt/research/compbio/krishnanlab/projects/chronic_inflammation/results/biogrid_noRWR/cluster_GOBP_enrichment/expanded_disease_genes_FDR=FALSE_cutoff=0.01_no_chron_inf_seed"
}

if(!dir.exists(outdir)) {dir.create(outdir)}

write.csv(all_results, file = paste0(outdir,
                                     "/",
                                     disease, 
                                     "_", 
                                     graph, 
                                     "_iterations=", 
                                     iterations, 
                                     "_FDR=",
                                     FDR,
                                     "_cutoff=",
                                     cutoff,
                                     "_GOBP_enrichment.csv"))

print("done")

    
    