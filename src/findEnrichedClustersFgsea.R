#' @param args[1] path to cluster output file
#' @param args[2] path to GenePlexus predictions for chronic inflammation genes
#' @param args[3] results path

# setup -------------------------------------------------------------------
library(tidyverse)
library(fgsea)
args <- commandArgs(TRUE)

# Make list of genes in each cluster (pathways argument)-------------------

# load cluster file
clusters = read.csv(args[1], row.names = 1)

# find for clusters with at least 5 genes
size5 = 
  clusters %>%
  group_by(Cluster) %>%
  tally(name = "nGenes") %>%
  filter(nGenes >= 5) %>%
  pull(Cluster)

# split into a list by cluster

gene_list = list()

for(clust in size5){
  
  gene_list[[clust]] = 
    clusters %>%
    filter(Cluster == clust) %>%
    pull(Gene)
    
} 

# Make vector of GenePlexus probabilities (stats argument)-----------------

genePlexus = read.delim(args[2])
genePlexus_info = basename(args[2])
genePlexus_info = gsub("--predictions.tsv", "", genePlexus_info)
genePLexus_name = gsub("--.*", "", genePlexus_info)
genePlexus_info = gsub(paste0(genePLexus_name, "--"), "", genePlexus_info)

probs = genePlexus$Probability
names(probs) = genePlexus$Entrez

# fgsea -------------------------------------------------------------------

res = 
  fgseaMultilevel(
    pathways = gene_list,
    stats = probs,
    minSize = 5,
    scoreType = "pos"
  )


# add clustering info -----------------------------------------------------
res_df = res[,1:6]
clust_info = 
  clusters %>%
  select(Disease,
         PredictionNetwork,
         PredictionFeatures,
         NegativesFrom,
         PredictionThreshold,
         Method,
         ClusterGraph,
         Resolution) %>%
  distinct()

res_df = cbind(res_df, clust_info)
res_df$ChronicInflammationSource = genePLexus_name
res_df$ChronicInfammationGenePlexus = genePlexus_info

# save --------------------------------------------------------------------
outdir = paste0(args[3], "/FGSEA")
if(!dir.exists(outdir)){dir.create(outdir)}

write.csv(res_df, file = paste0(outdir, 
                                "/", 
                                clust_info$Disease, 
                                "_",
                                clust_info$ClusterGraph,
                                "_", 
                                genePLexus_name,
                                "_fgsea_enrichment.csv"))


