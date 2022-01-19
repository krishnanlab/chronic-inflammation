# This file will submit up to 1000 qsubs at once
library("tidyverse")
#' @args[1] first index
#' @args[2] last index

args <- commandArgs(TRUE)
dirname <- "sbatches_permutedNetwork"
if(!dir.exists(dirname)){
  dir.create(dirname)
}

#first = as.numeric(args[1])
#print(paste0("first=", first))
#last = as.numeric(args[2])
#print(paste0("last=", last))

cluster_dirs = dir("/mnt/scratch/mckimale/chronic/biogrid_hotnet2_permutations/biogrid_hotnet2_permutations_noRWR/expanded_disease_genes_biogrid_iterations=2_FDR=FALSE_cutoff=0.01/clusters/",
full.names = FALSE)
cluster_dirs = cluster_dirs[2:length(cluster_dirs)]
cluster_dirs = gsub("biogrid_permutation_", "", cluster_dirs)
cluster_dirs = gsub("_ModularityVertexPartition_res=NA", "", cluster_dirs)
all = 1:1000
need = all[!all %in% as.numeric(cluster_dirs)]


#Edgelists remaining
edgelists=dir("/mnt/gs18/scratch/users/mckimale/chronic/biogrid_hotnet2_permutations/biogrid_hotnet2_permutations_noRWR/per_edgelists",
              full.names = FALSE)
edgelists=gsub("biogrid_permutation_","",edgelists)
edgelists=gsub("_edgelist_1","",edgelists)
all = 1:1000
need = all[!all %in% as.numeric(edgelists)]
for(i in need) {
  
  rjobsh <- paste0("permutedNetwork_", 
                   i,
                   ".rjob.sh"); cat(rjobsh, "\n")
  rjobConn <- file(paste0(dirname,"/",rjobsh))
  writeLines(c("#!/bin/sh -login",
               "#SBATCH --mem=30GB",
               paste0("#SBATCH --job-name=perm_", i),
               "#SBATCH --time=8:59:00",
               "#SBATCH --nodes=1",
               "#SBATCH --cpus-per-task=1",
               "#SBATCH --account=wang-krishnan",
               "",
               "cd ../bin/hotnet2/scripts",
               "",
               "module purge",
               "module load GCCcore/6.4.0 binutils/2.28 GNU/6.4.0-2.28 ncurses/6.1 CMake/3.11.1",
               ". /mnt/home/mckimale/anaconda2/etc/profile.d/conda.sh",
               "conda activate base",
               "which python",
               "",
               paste0("python2 permuteNetwork.py ",
                      "-e  ../../../data/biogrid_entrez_edgelist.txt ",
                      "-p  biogrid_permutation_",
                      i,
                      " -n  1 ",
                      "-o  /mnt/scratch/mckimale/chronic/biogrid_hotnet2_permutations/biogrid_hotnet2_permutations_noRWR ",
                      "-c  1")
               ),rjobConn)
  
  close(rjobConn)
  
  system(paste0("sbatch ", paste0(dirname,"/",rjobsh)))
  
  njobs <- system("squeue -u  mckimale | wc -l", intern=TRUE)
  njobs <- as.numeric(njobs)
  
  while(njobs > 1001) {
    Sys.sleep(360)
    njobs <- system("squeue -u mckimale | wc -l", intern=TRUE)
    njobs <- as.numeric(njobs)
  }
}




