#' @args[1] path to overlap scores files
#' @args[2] outdir
#' @args[3] real only T/F

args <- commandArgs(TRUE)
library(tidyverse)
library(parallel)
real_only_logical = as.logical(args[3])

#path = "/mnt/research/compbio/krishnanlab/projects/chronic_inflammation/results/biogrid_noRWR/expanded_disease_genes_biogrid_iterations=2_FDR=FALSE_cutoff=0.01/overlaps/scores"
path = args[1]
setwd(path)
files = list.files(path)

temp = mclapply(files, 
                read.csv, 
                row.names = 1,
                mc.cores= detectCores())

df = do.call(rbind, temp)
rm(temp)  

all_disease = unique(df$Disease)

ai = c("Celiac_Disease",
       "Crohn_Disease",
       "Irritable_Bowel_Syndrome",
       "Lupus_Erythematosus_Systemic",
       "Multiple_Sclerosis",
       "Psoriasis",
       "Rheumatoid_Arthritis",
       "Ulcerative_Colitis",
       "Vasculitis",
       "Diabetes_Mellitus_NonInsulinDependent"
)

neg =  c("Right_Hand_Phone",
         "North_UK",
         "Left_Hand_Phone",
         "Height",
         "Equal_Hand_Phone",
         "East_UK",
         "Driving_Speed",
         "Banana_Intake")

neg_controls = all_disease[grep(paste(neg, collapse="|"),
                                all_disease)]

fake_controls = all_disease[grep("Faketrait",
                                 all_disease)]


df$Control = ifelse(df$Disease %in% ai, 
                    "Autoimmune_disease", 
                    ifelse(df$Disease %in% neg_controls,
                           "Negative_trait",
                           ifelse(df$Disease %in% fake_controls,
                                  "Fake_trait",
                                  "Complex_disease")))

df$Permuted = ifelse(df$Graph == "biogrid", 
                     "real", 
                     "permuted")

permuted_df = df %>%
  filter(Permuted == "permuted")
colnames(permuted_df) = paste0(colnames(permuted_df), "_permuted")

mutatePermutedPval <-function(row,Group1, 
                              Disease, 
                              Enrichment, 
                              Permuted){
  
  print(paste0("current row = ", row))
  permuted_scores = permuted_df %>%
    filter(Group1_permuted == Group1[row], 
           Disease_permuted == Disease[row]) %>%
    pull(Enrichment_permuted)
  
  if(Permuted[row] == "real"){
    permuted_scores = c(permuted_scores, Enrichment[row]) 
  }
  
  num_permutation = length(permuted_scores)
  
  pval = sum(permuted_scores >= Enrichment[row])/num_permutation
  print(pval)
  return(pval)
}

if(real_only_logical == TRUE){
  df = filter(df, Permuted == "real") 
}

rowtoget=1:nrow(df)
print(paste0("nrows=", nrow(df)))
ans =  mclapply(rowtoget,
                Group1=df$Group1,
                Disease=df$Disease,
                Enrichment=df$Enrichment,
                Permuted=df$Permuted,
                mutatePermutedPval,
                mc.cores= detectCores())

df$PermutedPval = unlist(ans)

df = df %>%
  group_by(ChronicInf, 
           Disease, 
           Graph) %>%
  mutate(PermutedFDR = p.adjust(PermutedPval))

#outdir = "/mnt/research/compbio/krishnanlab/projects/chronic_inflammation/results/biogrid_noRWR/expanded_disease_genes_biogrid_iterations=2_FDR=FALSE_cutoff=0.01/overlaps"
outdir = args[2]

if(real_only_logical == TRUE){
  save(df, file = paste0(outdir, "/overlap_results_real_only.Rdata"))
  print("overlap_results_real_only.Rdata saved")
} else {
  save(df, file = paste0(outdir, "/overlap_results_real_and_permuted.Rdata"))
  print("overlap_results_real_and_permuted.Rdata saved")
}


    