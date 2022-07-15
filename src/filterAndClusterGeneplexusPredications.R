#' @param args[1] path to gene plexus predictions for real disease genes
#' @param args[2] prediction threshold
#' @param args[3] path to igraph object containing network for clustering
#' @param args[4] partition type
#' @param args[5] resolution parameter
#' @param args[6] results path
#' @param args[7] weighted T/F

# setup -------------------------------------------------------------------
library(tidyverse)
library(mccf1)
library(igraph)
# https://github.com/cole-trapnell-lab/leidenbase
library(leidenbase)
source("../src/chronic_inflammation_functions.R")
args <- commandArgs(TRUE)

# load disease gene predictions
print("load disease gene predictions")
pred = read.delim(args[1])

# get prediction file info
file_sep = str_split(basename(args[1]), pattern = "--")
file_sep = unlist(file_sep)
disease = file_sep[1]
pred_net = file_sep[2]
features = file_sep[3]
GSC = file_sep[4]
print(paste0(disease, "--", pred_net, "--", features))

# load network
g = loadRData(args[3])
clust_net = gsub("_igraph.Rdata", "", basename(args[3]))

#results_path = "/mnt/gs18/scratch/users/hickeys6/chronic_inflammation/GenePlexus_SLA_string"
results_path = args[6]
if(!dir.exists(results_path)) {dir.create(results_path)}

outdir = paste0(results_path,
                "/clusters")
if(!dir.exists(outdir)) {dir.create(outdir)}

outdir = paste0(outdir, "/predicted_with", pred_net, "--clustered_on_", clust_net)
if(!dir.exists(outdir)) {dir.create(outdir)}

# pick prediction threshold with mccf1 -------------------------------------
print("define prediction threshold")
if(args[2] == "mccf1"){
  response = 
    pred %>%
    filter(!Class.Label == "U") %>%
    mutate(Response = 
             case_when(Class.Label == "P" ~ 1,
                       Class.Label == "N" ~ 0)) %>%
    pull(Response)
  
  predictor = 
    pred %>%
    filter(!Class.Label == "U") %>%
    pull(Probability)
  
  mccf1_object = mccf1(response, predictor)
  mccf1_summary = summary(mccf1_object)
  threshold = mccf1_summary[1,2]
} else {
  threshold = as.numeric(args[2])
}

# filter pred -------------------------------------------------------------
print("filter disease genes")
disease_genes = 
  pred %>%
  filter(Probability >= threshold) %>%
  pull(Entrez)

# which are seed genes?
seeds = pred %>% filter(Class.Label == "P") %>% pull(Entrez)
# which are new genes?
new_genes = disease_genes[!disease_genes %in% seeds]

# nodes -------------------------------------------------------------------
keep_nodes = as.character(disease_genes)

# edges -------------------------------------------------------------------
# filter for edges incident on nodes from above
print("make disease subgraph")
sub.g = induced_subgraph(g, V(g)[name %in% keep_nodes])

# cluster -----------------------------------------------------------------
print("cluster disease genes")

if(as.logical(args[7]) == T){
  
  clusters = leiden_find_partition(sub.g,
                                   partition_type = args[4],
                                   #partition_type = "ModularityVertexPartition",
                                   resolution_parameter = as.numeric(args[5]),
                                   #resolution_parameter = .01,
                                   seed = 1,
                                   edge_weights = E(sub.g)$weight,
                                   num_iter = 100
                                   )
} else {
  
  clusters = leiden_find_partition(sub.g,
                                   partition_type = args[4],
                                   #partition_type = "ModularityVertexPartition",
                                   resolution_parameter = as.numeric(args[5]),
                                   #resolution_parameter = .01,
                                   seed = 1,
                                   num_iter = 100)
  
}

# make cluster df ---------------------------------------------------------
print("build cluster_df")
cluster_df = data.frame(Gene = as_ids(V(sub.g)), 
                        Disease = disease,
                        PredictionNetwork = pred_net,
                        PredictionFeatures = features,
                        NegativesFrom = GSC,
                        PredictionThreshold = args[2],
                        Cluster = paste0(disease, 
                                         "_",
                                         clusters$membership),
                        ClusterID = clusters$membership,
                        Method = args[4],
                        #Method = "ModularityVertexPartition",
                        ClusterGraph = clust_net)

if(args[4] == "ModularityVertexPartition"){
  cluster_df$Resolution = "NA"
  res = "NA"
} else {
  cluster_df$Resolution = as.numeric(args[5])
  res = args[5]
}

cluster_df$GeneType = 
  case_when(
    cluster_df$Gene %in% seeds ~ "seed",
    cluster_df$Gene %in% new_genes ~ "predicted")

# save cluster df --------------------------------------------------------------------
if(!dir.exists(outdir)) {dir.create(outdir)}

write.csv(cluster_df, 
          file = paste0(outdir,
                        "/", 
                        disease, 
                        "--threshold--",
                        args[2],
                        "--PredictionGraph--",
                        pred_net,
                        "--ClusterGraph--",
                        clust_net,
                        "_clusters.csv"))

print(paste0(outdir,
             "/", 
             disease, 
             "--threshold--",
             args[2],
             "--PredictionGraph--",
             pred_net,
             "--ClusterGraph--",
             clust_net,
             "_clusters.csv saved"))
