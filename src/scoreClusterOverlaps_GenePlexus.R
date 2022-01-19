#' @args[1] path to folder containing leiden cluster output files
#' @args[2] path to chronic inflammation prediction file
#' @args[3] path to output dir
#' @args[4] disease of interest
#' @args[5] chronic inflammation prediction threshold
#' @args[6] genes in network the disease genes were clustered on 
#' a one column csv file

library(tidyverse)
library(parallel)
args <- commandArgs(TRUE)
source("chronic_inflammation_functions.R")

output_dir = args[3]
if(!dir.exists(output_dir)){
  dir.create(output_dir)}

# chronic inflammation predictions -------------------------------------------
print("chronic inflammation clusters")
cioi = read.delim(args[2])
cioi_name = sub("--.*$", "", basename(args[2]))

# disease clusters --------------------------------------------------------
print("disease clusters")
disease_dir = args[1]
setwd(disease_dir)
# dir = "/mnt/research/compbio/krishnanlab/projects/chronic_inflammation/results/leiden_ModularityVertexPartition_b=0.5"
doi = args[4]

doi_grep = paste0("^", doi)
doi_fake_grep = paste0("Fake_", doi)

real_files = list.files(disease_dir, pattern = doi_grep)
fake_files = list.files(disease_dir, pattern = doi_fake_grep)
files = c(real_files, fake_files)

# read in the cluster files
cluster_list = lapply(files, read.csv, row.names = 1)
cluster_list = lapply(cluster_list, function(x) {
  y = 
    x %>% 
    dplyr::select(Gene, 
                  Disease, 
                  Cluster,
                  ClusterGraph); 
  y})

dis = do.call(rbind, cluster_list)
graph = unique(dis$ClusterGraph)

# find n genes per cluster
genes_per_cluster = dis %>%
  group_by(Cluster, Disease) %>%
  tally() %>%
  arrange(desc(n))

# filter out clusters with fewer than 5 genes
keep= genes_per_cluster %>%
  filter(n >= 5) %>%
  pull(Cluster)

if(length(keep) < 1){stop("no clusters with at least 5 genes")}

dis_filt = dis %>%
  filter(Cluster %in% keep)

all_diseases = unique(dis_filt$Disease)

# hypergeometric test -----------------------------------------------------
print("hypergeometric test")

network_genes = read.csv(args[6], row.names = 1)
network_genes = as.character(network_genes[,1])

ci_hyper = 
  cioi %>%
  filter(Probability >= as.numeric(args[5]),
         Entrez %in% network_genes) %>%
  mutate(Disease = cioi_name) %>%
  dplyr::select(Entrez, Disease)

colnames(ci_hyper) = c("Gene", "Cluster")

CI_genes  = 
  ci_hyper %>%
  pull(Gene)

scores = list()
shared_genes = list()

for(disease in all_diseases){
  
  print(disease)
  
  dis_hyper = 
    dis_filt %>%
    filter(Disease == disease) %>%
    dplyr::select(Gene, Cluster)
  
  dis_genes = 
    dis %>%
    filter(Disease == disease) %>%
    pull(Gene)

  bg = length(unique(c(CI_genes, dis_genes)))
  print(paste0("background = ", bg))
  
  overlaps = overlapSets(ci_hyper, dis_hyper, bg)
  print("overlapSets done")
  
  scores[[disease]] = overlaps[[2]]
  scores[[disease]]$Disease = sub("_[^_]+$", "", scores[[disease]]$Group2)
  scores[[disease]]$ChronicInf = sub("_[^_]+$", "", scores[[disease]]$Group1)
  row.names(scores[[disease]]) = NULL

  shared_genes[[disease]] = overlaps[[1]]
  row.names(shared_genes[[disease]]) = NULL
}

print("loop done")

score_dir = paste0(output_dir, "/scores")
if(!dir.exists(score_dir)){
  dir.create(score_dir)}

all_scores = do.call(rbind, scores)

all_scores = 
  all_scores %>%
  group_by(Disease) %>%
  mutate(FDR = p.adjust(Pval, "BH"))

all_shared_genes = do.call(rbind, shared_genes)

write.csv(all_scores, file = paste0(score_dir, 
                                "/", 
                                doi,
                                "_",
                                cioi_name, 
                                "_", 
                                graph, 
                                "_enrichment_scores.csv"))

print(paste0(score_dir, 
             "/", 
             doi, 
             "_",
             cioi_name, 
             "_", 
             graph, 
             "_enrichment_scores.csv saved"))

shared_dir = paste0(output_dir, "/shared_genes")
if(!dir.exists(shared_dir)){
  dir.create(shared_dir)}

write.csv(all_shared_genes, file = paste0(shared_dir, 
                                      "/", 
                                      doi,
                                      "_",
                                      cioi_name, 
                                      "_", 
                                      graph, 
                                      "_shared_genes.csv"))

print(paste0(shared_dir, 
             "/", 
             doi,
             "_",
             cioi_name, 
             "_", 
             graph, 
             "_shared_genes.csv saved"))

print("done")
