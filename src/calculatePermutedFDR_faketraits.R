#' @args[1] path to overlap scores files
#' @args[2] outdir
#' @args[3] real only T/F

args <- commandArgs(TRUE)
library(tidyverse)
library(parallel)
real_only_logical = as.logical(args[3])

path = args[1]
files = list.files(path, pattern = ".csv",full.names=T)

temp = mclapply(files, 
                read.csv, 
                row.names = 1,
                mc.cores= detectCores()-1)

df = do.call(rbind, temp)
rm(temp)  

df = 
  df %>%
  mutate(TraitType =
           case_when(
             grepl("Fake", Disease) == TRUE ~ "Fake",
             grepl("Fake", Disease) == FALSE ~ "Real"))

fake_df = 
  df %>%
  filter(TraitType == "Fake")

colnames(fake_df) = paste0(colnames(fake_df), "_fake")

fake_df$Disease_fake = gsub("Fake_", "", fake_df$Disease_fake)
fake_df$Disease_fake = gsub("_[^_]+$", "", fake_df$Disease_fake)

npermutations_df = 
  fake_df %>%
  group_by(Disease_fake) %>%
  tally %>%
  arrange(n)

colnames(npermutations_df) = c("Disease", "nFakeClusters")

mutatePermutedPval <-function(row,
                              Group1, 
                              Disease, 
                              Enrichment, 
                              TraitType){
  
  print(paste0("current row = ", row))
  permuted_scores = 
    fake_df %>%
    filter(Group1_fake == Group1[row], 
           Disease_fake == Disease[row]) %>%
    pull(Enrichment_fake)
  
  if(TraitType[row] == "Real"){
    permuted_scores = c(permuted_scores, Enrichment[row]) 
  }
  
  num_permutation = length(permuted_scores)
  
  pval = sum(permuted_scores >= Enrichment[row])/num_permutation
  print(pval)
  return(pval)
}

if(real_only_logical == TRUE){
  df = filter(df, TraitType == "Real") 
  
}

rowtoget=1:nrow(df)
print(paste0("nrows=", nrow(df)))
ans =  mclapply(rowtoget,
                Group1=df$Group1,
                Disease=df$Disease,
                Enrichment=df$Enrichment,
                TraitType=df$TraitType,
                mutatePermutedPval,
                mc.cores= detectCores())

df$PermutedPval = unlist(ans)

df = 
  df %>%
  group_by(ChronicInf, 
           Disease) %>%
  mutate(PermutedFDR = p.adjust(PermutedPval))

df = left_join(df, npermutations_df)

outdir = args[2]

if(real_only_logical == TRUE){
  save(df, file = paste0(outdir, "/overlap_results_real_only.Rdata"))
  print("overlap_results_real_only.Rdata saved")
} else {
  save(df, file = paste0(outdir, "/overlap_results_real_and_permuted.Rdata"))
  print("overlap_results_real_and_permuted.Rdata saved")
}
