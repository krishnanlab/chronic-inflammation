#' @args[1] path to overlap_results.Rdata
#' @args[2] cutoff
#' @args[3] output dir
#' @args[4] path to number_clustered.txt
#' @args[5] Path to gene cluster assignment files 

args <- commandArgs(TRUE)
source("/mnt/research/compbio/krishnanlab/projects/chronic_inflammation/src/chronic_inflammation_functions.R")
library(tidyverse)
library(parallel)

df = loadRData(args[1])
cutoff = as.numeric(args[2])
outdir = args[3]

# get rid of traits where there were clustering errors with the fake traits
# cd /mnt/research/compbio/krishnanlab/projects/chronic_inflammation/run/sbatches_clusterRandomGenes
# grep -hnr --with-filename "number clustered" *.out > number_clustered.txt

nclust = read.delim(args[4], header = F, sep = " ")
nclust$nFakesClustered = gsub("number clustered = ", "", nclust$V2)
nclust$nFakesClustered = as.numeric(nclust$nFakesClustered)

nclust$Disease = gsub("_BioGrid.*", "", nclust$V1)

remove = 
  nclust %>%
  filter(nFakesClustered < 5000) %>%
  pull(Disease)

# fix it for alex
final_alex = df %>%
  filter(PermutedFDR <= cutoff,
         !Disease %in% remove)

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
  "phyperFDR",
  "PermutedPval",
  "PermutedFDR")

final_alex = 
  final_alex %>%
  ungroup() %>%
  select(Chronic_Inflamation_Cluster,
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

write.csv(final_alex, file = paste0(outdir, "/final_for_alex.csv"))

# make df with gene assignments from real and fake clusters
cluster_path = args[5]
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

write.csv(df_filt, file = paste0(outdir, "/relevant_gene_cluster_assigments.csv"))
          
    
          
          











