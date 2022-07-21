#' @args[1] path to overlap_results.Rdata
#' @args[2] cutoff
#' @args[3] output dir
# @args[4] path to number_clustered.txt
#' @args[4] Path to gene cluster assignment files 

args <- commandArgs(TRUE)
source("../src/chronic_inflammation_functions.R")
library(tidyverse)
library(parallel)

results_file = args[1]  
df = loadRData(results_file) 

results_file = basename(results_file)

ci_gene_set = gsub("_thresh.*", "", results_file)
clusterGraph = gsub("^.*clusteredOn_", "", results_file)
clusterGraph = gsub("_overlap_results.Rdata", "", clusterGraph)
predNet = gsub("^.*predicted_with_", "", results_file)
predNet = gsub("_clusteredOn.*", "", predNet)

cutoff = as.numeric(args[2])
outdir = args[3]

# fix it for alex
final_alex = df %>%
  filter(PermutedFDR <= cutoff)
         #!Disease %in% remove)

colnames(final_alex) = c(
  "Chronic_Inflamation_Cluster",
  "n_Chronic_Inflamation_Cluster_genes",
  "Disease_Cluster",
  "n_Disease_Cluster_genes",
  "n_overlap",
  "phyperPval",
  "Enrichment",
  "Disease",
  "ChronicInf",
  "ChronicInfThreshold",
  "ClusterGraph",
  "PredictionNetwork",
  "phyperFDR",
  "PermutedPval",
  "PermutedFDR",
  "nFakeClusters")

final_alex = 
  final_alex %>%
  ungroup() %>%
  select(PredictionNetwork,
         ClusterGraph,
         Chronic_Inflamation_Cluster,
         Disease_Cluster,
         n_overlap,
         n_Disease_Cluster_genes,
         n_Chronic_Inflamation_Cluster_genes,
         Enrichment,
         Disease,
         phyperPval,
         phyperFDR,
         PermutedPval,
         PermutedFDR)

write.csv(final_alex, 
          file = paste0(outdir, 
                        "/", 
                        ci_gene_set,
                        "--predictedWith--",
                        predNet,
                        "--clusteredOn--", 
                        clusterGraph, 
                        "_final_for_alex.csv"))

# make df with gene assignments from real and fake clusters
cluster_path = args[4]
files = list.files(cluster_path, full.names = T, pattern = ".csv")

doi = unique(final_alex$Disease)
files_oi = grep(paste(doi,collapse="|"), files, value = TRUE)

temp = mclapply(files_oi,
                read.csv, 
                row.names = 1,
                mc.cores = detectCores()-1)

temp = 
  mclapply(temp, function(x) {
    
    y = 
      x %>% 
      dplyr::select(Gene, 
                    Disease, 
                    Cluster,
                    ClusterGraph); 
    y},
    
    mc.cores = detectCores()-1)

df = do.call(rbind, temp)

cluster_size = 
  df %>%
  group_by(Cluster) %>%
  tally()

greater5 = 
  cluster_size %>%
  filter(n >= 5) %>%
  pull(Cluster)

df_filt =
  df %>%
  filter(Cluster %in% greater5)

write.csv(df_filt, 
          file = paste0(outdir, 
                        "/",
                        ci_gene_set,
                        "--predictedWith--",
                        predNet,
                        "--clusteredOn--", 
                        clusterGraph, 
                        "_relevant_gene_cluster_assigments.csv"))
          
    
          
          











