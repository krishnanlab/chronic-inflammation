#' @param args[1] path to original gene file (entrez ids)
#' @param args[2] path to list of node neighbors from edgelist
#' @param args[3] path to edgelist
#' @param args[4] number of hypergeometric iterations
#' @param args[5] pval cutoff
#' @param args[6] FDR = TRUE, pval = FALSE
#' @param args[7] partition_type
#' @param args[8] resolution parameter
#' @param args[9] results path

# setup -------------------------------------------------------------------
library(tidyverse)
library(igraph)
# https://github.com/cole-trapnell-lab/leidenbase
library(leidenbase)
source("/mnt/research/compbio/krishnanlab/projects/chronic_inflammation/src/chronic_inflammation_functions.R")


args <- commandArgs(TRUE)

n_iterations = as.numeric(args[4])
cutoff = as.numeric(args[5])
FDR_logical = args[6]
netname = gsub("_neighbor_list.Rdata", "", basename(args[2]))

#results_path = "/mnt/research/compbio/krishnanlab/projects/chronic_inflammation/results/biogrid_noRWR"
results_path = args[9]

outdir = paste0(results_path,
                "/expanded_disease_genes_biogrid_iterations=",
                n_iterations,
                "_FDR=",
                FDR_logical,
                "_cutoff=", 
                cutoff)

if(!dir.exists(outdir)) {dir.create(outdir)}

# expand genes list --------------------------------------------------------
# load original genes
og = read.delim(args[1])
disease = gsub(".txt", "", basename(args[1]))
# load list of first degree neighbors
nh_list = loadRData(args[2])

# perform hypergeometric test between original genes and node neighborhoods
# keep nodes with significant overlap
# iterate
new_list = expandGeneList(og, 
                          nh_list, 
                          cutoff = cutoff,
                          FDR_TF = FDR_logical,
                          iterations = n_iterations)

new_nodes = new_list[[1]]
new_nodes$expandFDR = FDR_logical
new_nodes$expandCutoff = cutoff

n_iterations = new_list[[2]]
new_nodes$Iterations = n_iterations

# save the expanded gene df to a folder named after
# the network
dg_dir = paste0(outdir, "/disease_genes")
if(!dir.exists(dg_dir)) {dir.create(dg_dir)}

write.table(new_nodes, 
            file = paste0(dg_dir, 
                          "/", 
                          disease, 
                          ".txt"),
            sep = "\t",
            quote = F,
            row.names = F
)

# nodes -------------------------------------------------------------------
keep_nodes = as.character(new_nodes[,1])
# edges -------------------------------------------------------------------
edgelist = read.delim(args[3], header = F)
edgelist$V1 = as.character(edgelist$V1)
edgelist$V2 = as.character(edgelist$V2)

g = graph_from_edgelist(as.matrix(edgelist), directed = FALSE)

# filter for edges incident on nodes from above
sub.g = induced_subgraph(g, keep_nodes)

# cluster -----------------------------------------------------------------
clusters = leiden_find_partition(sub.g,
                                 partition_type = args[7],
                                 resolution_parameter = as.numeric(args[8]),
                                 seed = 1)

# make cluster df ---------------------------------------------------------
cluster_df = data.frame(Gene = as_ids(V(sub.g)), 
                        Disease = disease,
                        Cluster = paste0(disease, 
                                         "_",
                                         clusters$membership),
                        ClusterID = clusters$membership,
                        Method = args[7],
                        #Method = "ModularityVertexPartition",
                        Graph = netname,
                        Iterations = n_iterations,
                        expandFDR = FDR_logical,
                        expandCutoff = cutoff)

if(args[7] == "ModularityVertexPartition"){
  cluster_df$Resolution = "NA"
  res = "NA"
} else {
  cluster_df$Resolution = as.numeric(args[8])
  res = args[8]
}

# save disease graph --------------------------------------------------------------------
graph_dir =  paste0(outdir, "/disease_graphs")
if(!dir.exists(graph_dir)) {dir.create(graph_dir)}
save(sub.g, file = paste0(graph_dir, 
                             "/", disease, 
                             "_igraph.Rdata"))

# save cluster df --------------------------------------------------------------------
cluster_dir = paste0(outdir, "/clusters")
if(!dir.exists(cluster_dir)) {dir.create(cluster_dir)}

write.csv(cluster_df, 
          file = paste0(cluster_dir,
                        "/", 
                        disease, 
                        "_clusters.csv"))

# for testing best parameters --------------------------------------------------------
# per disease

parameter_check_dir = paste0(results_path, "/parameter_checks")
if (!dir.exists(parameter_check_dir)) {dir.create(parameter_check_dir)}

min_size = 5

keep_clusters = cluster_df %>%
  group_by(Cluster) %>%
  tally() %>%
  filter(n >=min_size)

print(paste0("number of clusters with >= 5 genes = ", nrow(keep_clusters)))

total_genes = nrow(cluster_df)
kept_genes = sum(keep_clusters$n)
frac = kept_genes/total_genes

n_clusters = nrow(keep_clusters)

per_disease = data.frame(Disease = disease,
                         Graph = netname,
                         Method = args[7],
                         Resolution = res,
                         Iterations = n_iterations,
                         expandFDR = FDR_logical,
                         expandCutoff = cutoff,
                         Modularity = clusters$modularity,
                         MinSize = min_size,
                         totalGenes = total_genes,
                         clusteredGenes = kept_genes,
                         FractionTotalGenes = frac,
                         nClusters = n_clusters)

write.csv(per_disease, 
          file = paste0(
            parameter_check_dir,
            "/",
            disease,
            "_n_iterations=",
            n_iterations,
            ".csv"),
          row.names = F)

# per cluster
# save per cluster output
leiden_dir = paste0(results_path, 
                    "/leiden_alg_per_cluster_output")

if (!dir.exists(leiden_dir)) {dir.create(leiden_dir)}

nGenes = cluster_df %>%
  group_by(Cluster, ClusterID) %>%
  tally() %>%
  arrange(as.numeric(ClusterID)) 

per_cluster = data.frame(Disease = disease,
                         Graph = netname,
                         Method = args[7],
                         Resolution = res,
                         Iterations = n_iterations,
                         expandFDR = FDR_logical,
                         expandCutoff = cutoff,
                         Cluster = nGenes$Cluster,
                         nGenes = nGenes$n,
                         edgeWeightWITHINcommunity = 
                           clusters$edge_weight_within_community,
                         edgeWeightFROMcommunity = 
                           clusters$edge_weight_from_community,
                         edgeWeightTOcommunity = 
                           clusters$edge_weight_to_community,
                         totalEdgeWeight = clusters$total_edge_weight 
                         )

write.csv(per_cluster, 
          file = paste0(leiden_dir, 
                        "/", 
                        disease, 
                        "_iterations=",
                        n_iterations,
                        ".csv"), 
          row.names = F)

print("done")
