#' @param args[1] path to cluster file
#' @param args[2] background genes from network, 1 col df, csv
#' @param args[3] outdir

# setup -------------------------------------------------------------------
library(topGO)
library(org.Hs.eg.db)
library(tidyverse)
source("../src/chronic_inflammation_functions.R")

args <- commandArgs(TRUE)
cluster_df = read.csv(args[1], row.names = 1)
network_genes = read.csv(args[2], row.names = 1)
network_genes = as.character(network_genes[,1])

meta = cluster_df %>%
  dplyr::select(ChronicInflammationSource, 
                Method,
                PredictionNetwork,
                ClusterGraph,
                Resolution) %>%
  distinct()

keep_clusters = cluster_df %>%
  group_by(Cluster) %>%
  tally() %>%
  filter(n >=5) %>%
  pull(Cluster)

clusters = unique(cluster_df$Cluster)
clusters = clusters[clusters %in% keep_clusters]
all_results = {}
if(length(keep_clusters) < 1){
  print("no clusters with >= 5 genes")
} else {
  
  for(clust in clusters){
    
    print(clust)
    
    clust_genes = cluster_df %>%
      filter(Cluster == clust) %>%
      pull(Gene)
    
    background = network_genes
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
  
  outdir = args[3]
  
  if(!dir.exists(outdir)) {dir.create(outdir)}
  
  write.csv(all_results, file = paste0(outdir,
                                       "/",
                                       meta$ChronicInflammationSource,
                                       "--predictedWith--",
                                       meta$PredictionNetwork,
                                       "--clusteredOn--",
                                       meta$ClusterGraph,
                                       "--",
                                       meta$Method,
                                       "--resolution--",
                                       meta$Resolution,
                                       "_GOBP_enrichment.csv"))
  
  print("done")
  
}