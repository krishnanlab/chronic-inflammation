#' @param args[1] path to file with all random disease gene sets
#' @param args[2] disease 
#' @param args[3] path to igraph object containing network for clustering
#' @param args[4] partition type
#' @param args[5] resolution parameter
#' @param args[6] results path

# setup -------------------------------------------------------------------
library(tidyverse)
library(igraph)
# https://github.com/cole-trapnell-lab/leidenbase
library(leidenbase)
source("chronic_inflammation_functions.R")
args <- commandArgs(TRUE)

# load random
print("load disease gene predictions")
all_random = read.delim(args[1])

doi = args[2]
doi_grep = paste0("Fake_",doi)

g = loadRData(args[3])
netname = gsub("_igraph.Rdata", "", basename(args[3]))

outdir = args[6]
if(!dir.exists(outdir)) {dir.create(outdir)}

doi_random = 
  all_random %>%
  filter(grepl(doi_grep, Disease))

random_list = 
  doi_random %>%
  group_split(Disease) %>%
  setNames(unique(doi_random$Disease))

number_clustered = 0

all_clusters = {}
for(df in random_list){
  
  disease = unique(df$Disease)
  print(disease)
  keep_nodes = as.character(df$Gene)
  
  # filter for edges incident on nodes from above
  print("make disease subgraph")
  sub.g = induced_subgraph(g, V(g)[name %in% keep_nodes])
  
  print("cluster disease genes")
  clusters = try(leiden_find_partition(sub.g,
                                   partition_type = args[4],
                                   #partition_type = "ModularityVertexPartition",
                                   resolution_parameter = as.numeric(args[5]),
                                   #resolution_parameter = .01,
                                   seed = 1))
  
  if(class(clusters) == "try-error"){print("leiden error")} else {number_clustered = number_clustered + 1}
  if(class(clusters) == "try-error"){next}
  
  print("build cluster_df")
  cluster_df = data.frame(Gene = as_ids(V(sub.g)), 
                          Disease = disease,
                          Cluster = paste0(disease, 
                                           "_",
                                           clusters$membership),
                          ClusterID = clusters$membership,
                          Method = args[4],
                          #Method = "ModularityVertexPartition",
                          ClusterGraph = netname)
  
  if(args[4] == "ModularityVertexPartition"){
    cluster_df$Resolution = "NA"
    res = "NA"
  } else {
    cluster_df$Resolution = as.numeric(args[5])
    res = args[5]
  }
  
  all_clusters = rbind(all_clusters, cluster_df)
}

write.csv(all_clusters, 
          file = paste0(outdir,
                        "/Fake_", 
                        doi, 
                        "--ClusterGraph--",
                        netname,
                        "_clusters.csv"))

print(paste0(outdir,
             "/Fake_", 
             doi, 
             "--ClusterGraph--",
             netname,
             "_clusters.csv saved"))

print(paste("number clustered =", number_clustered))

numclusdir="../results/num_random_clustered"
if(!dir.exists(numclusdir)) {dir.create(numclusdir)}
write_tsv(tibble(Disease=doi,nFakesClustered=number_clustered),paste0(numclusdir,doi,".txt"))