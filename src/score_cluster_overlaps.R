#' @args[1] path to folder containing leiden cluster output files
#' @args[2] path to real chronic inflammation clusters of interest
#'          c("chronic_inflammation_gene_shot_pubs_greater10.txt",
#'            "chronic_inflammation_gene_shot.txt",
#'            "chronic_inflammation_go.txt")
#' @args[3] path to output dir
#' @args[4] number of background genes

 
library(tidyverse)
args <- commandArgs(TRUE)
source("/mnt/research/compbio/krishnanlab/projects/chronic_inflammation/src/chronic_inflammation_functions.R")

output_dir = args[3]
if(!dir.exists(output_dir)){
  dir.create(output_dir)}

# chronic inflammation clusters -------------------------------------------
print("chronic inflammation clusters")
cioi = read.csv(args[2], row.names = 1)
cioi_name = sub("_clusters.csv", "", basename(args[2]))

# find n genes per cluster
genes_per_cluster = cioi %>%
  group_by(Cluster, Disease) %>%
  tally() %>%
  arrange(desc(n))

# filter out clusters with fewer than 5 genes
keep= genes_per_cluster %>%
  filter(n >= 5) %>%
  pull(Cluster)

cioi_filt = cioi %>%
  filter(Cluster %in% keep)

# disease clusters --------------------------------------------------------
print("disease clusters")
disease_dir = args[1]
setwd(disease_dir)
# dir = "/mnt/research/compbio/krishnanlab/projects/chronic_inflammation/results/leiden_ModularityVertexPartition_b=0.5"
files = list.files(disease_dir, pattern = ".csv")

# remove chronic inflammation clusters we don't care about
ci_files = c("chronic_inflammation_gene_shot_pubs_greater10_clusters.csv",
             "chronic_inflammation_gene_shot_clusters.csv",
             "chronic_inflammation_go_clusters.csv") 

files = files[!files %in% ci_files]

# read in the cluster files
cluster_list = lapply(files, read.csv, row.names = 1)

# bind all diseases into one df
dis = do.call(rbind, cluster_list)

# find n genes per cluster
genes_per_cluster = dis %>%
  group_by(Cluster, Disease) %>%
  tally() %>%
  arrange(desc(n))

# filter out clusters with fewer than 5 genes
keep= genes_per_cluster %>%
  filter(n >= 5) %>%
  pull(Cluster)

dis_filt = dis %>%
  filter(Cluster %in% keep)

# hypergeometric test -----------------------------------------------------
print("hypergeometric test")
ci_hyper = cioi_filt %>%
  dplyr::select(Gene, Cluster)

dis_hyper = dis_filt %>%
  dplyr::select(Gene, Cluster)

bg = as.numeric(args[4])

overlaps = overlapSets(ci_hyper, dis_hyper, bg)
print("overlapSets done")

scores = overlaps[[2]]
scores$Disease = sub("_[^_]+$", "", scores$Group2)
scores$ChronicInf = sub("_[^_]+$", "", scores$Group1)

scores = scores %>%
  group_by(Disease) %>%
  mutate(FDR = p.adjust(Pval, "BH"))

shared_genes = overlaps[[1]]

graph = unique(dis$Graph)
scores$Graph = graph
shared_genes$Graph = graph

score_dir = paste0(output_dir, "/scores")
if(!dir.exists(score_dir)){
  dir.create(score_dir)}

write.csv(scores, file = paste0(score_dir, 
                               "/", 
                               cioi_name, 
                               "_", 
                               graph, 
                               "_enrichment_scores.csv"))

print(paste0(score_dir, 
             "/", 
             cioi_name, 
             "_", 
             graph, 
             "_enrichment_scores.csv saved"))

shared_dir = paste0(output_dir, "/shared_genes")
if(!dir.exists(shared_dir)){
  dir.create(shared_dir)}

write.csv(shared_genes, file = paste0(shared_dir, 
                                "/", 
                                cioi_name, 
                                "_", 
                                graph, 
                                "_shared_genes.csv"))

print(paste0(shared_dir, 
             "/", 
             cioi_name, 
             "_", 
             graph, 
             "_shared_genes.csv saved"))

print("done")
