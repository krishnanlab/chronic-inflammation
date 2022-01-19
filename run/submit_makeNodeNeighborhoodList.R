# This file will submit up to 1000 qsubs at once
library("tidyverse")
#' @args[1] first index
#' @args[2] last index
#args <- commandArgs(TRUE)
#first = as.numeric(args[1])
#print(paste0("first=", first))
#last = as.numeric(args[2])
#print(paste0("last=", last))

first=1
last=1000
# path(s) folder containing original edge list
# and hotnet2 generated influence matrix
path = "../../biogrid_hotnet2_permutations/biogrid_hotnet2_permutations_noRWR/per_edgelists"
files = list.files(path, pattern = "edgelist_1")

dirname <- "sbatches_makeNodeNeighborhoodList"
if(!dir.exists(dirname)){
  dir.create(dirname)
}

for(i in c(first:last)) {
  
  rjobsh <- paste0("list_biogrid_edgelist_", 
                   i,
                   "_neighbors",
                   ".rjob.sh"); cat(rjobsh, "\n")
  rjobConn <- file(paste0(dirname,"/",rjobsh))
  writeLines(c("#!/bin/sh -login",
               "#SBATCH --mem=30GB",
               paste0("#SBATCH --job-name=list", i),
               "#SBATCH --time=2:00:00",
               "#SBATCH --nodes=1",
               "#SBATCH --cpus-per-task=1",
               "#SBATCH --account=wang-krishnan",
               "",
               "cd ../src",
               "",
               "module load GCC/8.3.0  OpenMPI/3.1.4  R/4.0.2",
               "",
               paste0("Rscript makeNodeNeighborhoodList.R ",
                      path,
                      "/biogrid_permutation_",
                      i,
                      "_edgelist_1",
                      " ",
                      FALSE,
                      " ",
                      TRUE,
                      " ",
                      "biogrid_permutation",
                      i,
                      " ",
                      "../../biogrid_hotnet2_permutations/biogrid_hotnet2_permutations_noRWR/node_neighborhood_lists"
                      )),rjobConn)
  
  close(rjobConn)
  
  system(paste0("sbatch ", paste0(dirname,"/",rjobsh)))
  
  njobs <- system("squeue -u  mckimale | wc -l", intern=TRUE)
  njobs <- as.numeric(njobs)
  
  while(njobs > 1000) {
    Sys.sleep(360)
    njobs <- system("squeue -u mckimale | wc -l", intern=TRUE)
    njobs <- as.numeric(njobs)
  }
}
