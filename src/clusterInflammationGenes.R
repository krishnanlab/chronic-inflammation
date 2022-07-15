#' @param args[1] path to inflammation genes
#' @param args[2] path to igraph object containing network for clustering
#' @param args[3] partition type
#' @param args[4] resolution parameter
#' @param args[5] results path
#' @param args[6] weighted T/F

# setup -------------------------------------------------------------------
library(tidyverse)
library(igraph)
# https://github.com/cole-trapnell-lab/leidenbase
library(leidenbase)
source("../src/chronic_inflammation_functions.R")
args <- commandArgs(TRUE)

# load inflammation genes
print("load inflammation genes")
inflam_genes = read.delim(args[1])
inflam_name = gsub(".txt", "", basename(args[1]))

# load network
g = loadRData(args[2])
netname = gsub("_igraph.Rdata", "", basename(args[2]))

# results path
results_path = args[5]
if(!dir.exists(results_path)) {dir.create(results_path)}

outdir = paste0(results_path, "/clusters")
if(!dir.exists(outdir)) {dir.create(outdir)}

# nodes -------------------------------------------------------------------
keep_nodes = as.character(inflam_genes$Gene)

# edges -------------------------------------------------------------------
# filter for edges incident on nodes from above
print("make disease subgraph")
sub.g = induced_subgraph(g, V(g)[name %in% keep_nodes])

# cluster -----------------------------------------------------------------
print("cluster disease genes")

if(as.logical(args[6]) == T){
  
  clusters = leiden_find_partition(sub.g,
                                   partition_type = args[3],
                                   #partition_type = "ModularityVertexPartition",
                                   resolution_parameter = as.numeric(args[4]),
                                   #resolution_parameter = .01,
                                   seed = 1,
                                   edge_weights = E(sub.g)$weight)
} else {
  
  clusters = leiden_find_partition(sub.g,
                                   partition_type = args[3],
                                   #partition_type = "ModularityVertexPartition",
                                   resolution_parameter = as.numeric(args[4]),
                                   #resolution_parameter = .01,
                                   seed = 1)
  
}

# make cluster dfs ---------------------------------------------------------
print("build cluster_df")
cluster_df = data.frame(Gene = as_ids(V(sub.g)), 
                        ChronicInflammationSource = inflam_name,
                        Cluster = paste0(inflam_name, 
                                         "_",
                                         clusters$membership),
                        ClusterID = clusters$membership,
                        #Method = args[3],
                        Method = "ModularityVertexPartition",
                        ClusterGraph = netname)

if(args[3] == "ModularityVertexPartition"){
  cluster_df$Resolution = "NA"
  res = "NA"
} else {
  cluster_df$Resolution = as.numeric(args[4])
  res = args[4]
}
 # save
write.csv(cluster_df, 
          file = paste0(outdir,
                        "/", 
                        inflam_name, 
                        "--ClusterGraph--",
                        netname,
                        "_clusters.csv"))

# save cluster dfs ---------------------------------------------------------
# save in separate files for geneplexus
# only save clusters with at least 5 genes

clust5 = 
  cluster_df %>%
  group_by(Cluster) %>%
  tally() %>%
  filter(n >= 5) %>%
  pull(Cluster)

cluster_list = 
  cluster_df %>%
  filter(Cluster %in% clust5) %>%
  group_by(Cluster) %>%
  group_split()

names(cluster_list) = unlist(lapply(cluster_list, function(x) {unique(x$Cluster)}))

for(i in 1:length(cluster_list)){
  
  for_gp = 
    cluster_list[[i]] %>%
    select(Gene, Cluster)
  
  write.table(for_gp,
              file = paste0(outdir,
                            "/", 
                            names(cluster_list)[i], 
                            "--ClusterGraph--",
                            netname,
                            ".txt"),
              sep = "\t",
              quote  = F,
              row.names = F)
  
}

