#' @param args[1] path to folder containing original gene files (entrez ids)
#' @param args[2] path to list of node neighbors from edgelist
#' @param args[3] path to edgelist
#' @param args[4] number of hypergeometric iterations
#' @param args[5] pval cutoff
#' @param args[6] FDR = TRUE, pval = FALSE
#' @param args[7] partition_type
#' @param args[8] resolution parameter
#' @param args[9] results path

# setup
library(tidyverse)
library(igraph)
# https://github.com/cole-trapnell-lab/leidenbase
library(leidenbase)
source("chronic_inflammation_functions.R")
args <- commandArgs(TRUE)

files = list.files(args[1], full.names = TRUE)
og_list = lapply(files, read.delim)
names(og_list) = gsub(".txt", "", basename(files))

nh_list = loadRData(args[2])
# edges 
print("earfased")

#read in real with tab, read in permuted with " "
edgelist = read.delim(args[3], header = F,sep="\t")
#If only 1 column, then its a permutation so reload with correct separater
if(ncol(edgelist)==1){
  edgelist=read.delim(args[3], header = F,sep=" ")
}
edgelist$V3 = NULL
edgelist$V1 = as.character(edgelist$V1)
edgelist$V2 = as.character(edgelist$V2)


g = graph_from_edgelist(as.matrix(edgelist), directed = FALSE)

netname = gsub("_edgelist_1", "", basename(args[3]))
n_iterations = as.numeric(args[4])
cutoff = as.numeric(args[5])
FDR_logical = args[6]


if(args[7] == "ModularityVertexPartition"){
  res = "NA"
} else {
  res = args[8]
}

# prep outdirs
results_path = args[9]
outdir = paste0(results_path,
                "/expanded_disease_genes_biogrid_iterations=",
                n_iterations,
                "_FDR=",
                FDR_logical,
                "_cutoff=", 
                cutoff)
print(outdir)
if(!dir.exists(outdir)) {dir.create(outdir)}

dg_dir = paste0(outdir, "/disease_genes")
if(!dir.exists(dg_dir)) {dir.create(dg_dir)}

netname_dg_dir =  paste0(dg_dir, "/", netname)
if(!dir.exists(netname_dg_dir)) {dir.create(netname_dg_dir)}

cluster_dir = paste0(outdir, "/clusters")
if(!dir.exists(cluster_dir)) {dir.create(cluster_dir)}

netname_clust_dir =  paste0(cluster_dir, 
                            "/", 
                            netname,
                            "_", args[7], 
                            "_res=", 
                            res)

if(!dir.exists(netname_clust_dir)) {dir.create(netname_clust_dir)}

for(i in 1:length(og_list)){
  
  og = og_list[[i]]
  disease = names(og_list)[i]
  print(paste0("start ", disease))
  # perform hypergeometric test between original genes and node neighborhoods
  # keep nodes with significant overlap
  # iterate
  new_list = expandGeneList(og, 
                            nh_list, 
                            cutoff = as.numeric(args[5]),
                            FDR_TF = args[6],
                            iterations = as.numeric(args[4]))
  
  new_nodes = new_list[[1]]
  new_nodes$expandFDR = FDR_logical
  new_nodes$expandCutoff = cutoff
  
  n_iterations = new_list[[2]]
  new_nodes$Iterations = n_iterations
  
  write.table(new_nodes, 
              file = paste0(netname_dg_dir, 
                            "/", 
                            disease, 
                            ".txt"),
              sep = "\t",
              quote = F,
              row.names = F
  )
  
  # nodes
  keep_nodes = as.character(new_nodes[,1])
  sub.g = induced_subgraph(g, keep_nodes)
  
  #cluster
  clusters = leiden_find_partition(sub.g,
                                   partition_type = args[7],
                                   resolution_parameter = as.numeric(args[8]),
                                   seed = 1)
  # make cluster df 
  cluster_df = data.frame(Gene = as_ids(V(sub.g)), 
                          Disease = disease,
                          Cluster = paste0(disease, 
                                           "_",
                                           clusters$membership),
                          ClusterID = clusters$membership,
                          Method = args[7],
                          Resolution = res,
                          Graph = netname,
                          Iterations = n_iterations,
                          expandFDR = FDR_logical,
                          expandCutoff = cutoff)
  
  #save
  write.csv(cluster_df, file = paste0(netname_clust_dir, 
                                      "/", 
                                      disease, 
                                      "_clusters.csv"))
  print(paste0(disease, " done"))
}

print(paste0(netname, " all done!"))
