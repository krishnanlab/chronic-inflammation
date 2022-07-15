#' @param args[1] saverunner input dir path 
#' @param args[2] path to interactome edgelist
#' @param args[3] path to "final for alex" file with significant clusters
#' @param args[4] path to gene cluster assigments

args <- commandArgs(TRUE)
library(tidyverse)

path = args[1]
# interactome needs to be string exp
int = read.delim(args[2], header = F)

int = 
  int %>%
  select(V1, V2) %>%
  rename(Gene_A_Entrez_ID = V1,
         GeneB_Entrez_ID = V2)

# needs to match the saverunner example
colnames(int)[1] = "Gene_A_Entrez ID"

write.table(int, 
            file = paste0(path, "/interactome.txt"),
            row.names = F,
            quote = F,
            sep = "\t")

# both the drug and disease files will be CI enriched cluster genes
genes = read.csv(args[4])
sig = read.csv(args[3])

load("/mnt/research/compbio/krishnanlab/projects/chronic_inflammation/results/GenePlexus_parameter_checks/Geneplexus_summary_new_networks.Rdata")

cv_df =
  summary_df %>%
  filter(Network == "ConsensusPathDB") %>%
  dplyr::select(Disease, 
                nSeeds,
                CVscore)

keep_disease = 
  cv_df %>%
  filter(CVscore >= 1,
         nSeeds >= 15,
         !Disease == "Colitis") %>%
  pull(Disease)

keep_cluster = 
  sig %>%
  filter(Disease %in% keep_disease) %>%
  pull(Disease_Cluster) %>%
  unique()

Disease = 
  genes %>% 
  filter(Cluster %in% keep_cluster) %>%
  select(Disease, Gene) %>%
  rename(disease = Disease,
         GeneID = Gene)

write.table(Disease, 
            file = paste0(path, "/disgenes.txt"),
            row.names = F,
            quote = F,
            sep = "\t")

Drugs = 
  Disease %>%
  rename(Drug = disease)

write.table(Drugs, 
            file = paste0(path, "/drugs.txt"),
            row.names = F,
            quote = F,
            sep = "\t")
