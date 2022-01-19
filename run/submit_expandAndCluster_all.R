# This file will submit up to 1000 qsubs at once
#' @param args[1] first perm
#' @param args[2] last perm 
#' @param args[3] n_iterations
library("tidyverse")
args <- commandArgs(TRUE)

# path to folder with original disease gene list
#for args[1]
dg_path = "../results/disease_gene_files"
setwd(dg_path)
dg_files = list.files(pattern = ".txt")
# path(s) to list of node neighbors from edgelist
# for args[2]
neighbors_path = "../../biogrid_hotnet2_permutations/biogrid_hotnet2_permutations_noRWR/node_neighborhood_lists"
# path to permuted edgelists
# for args[3]
edgelist_path = "../../biogrid_hotnet2_permutations/biogrid_hotnet2_permutations_noRWR/per_edgelists"
# results path for args[9]
results_path = "../../biogrid_hotnet2_permutations/biogrid_hotnet2_permutations_noRWR"
# number of hypergeo iterations
n_iter = as.numeric(args[3])
# set pval cutoff
cutoff = 0.01
# use FDR?
FDR_logical = FALSE
# leiden alg partitian type
partition_type = "ModularityVertexPartition"
# resolution parameter
res = 0.1

setwd("../../run")

dirname <- "sbatches_expandAndCluster"
if(!dir.exists(dirname)){
  dir.create(dirname)
}




first_perm = as.numeric(args[1])
last_perm = as.numeric(args[2])
for(i in first_perm:last_perm) {
    print(getwd())
    
    rjobsh <- paste0("all_dsieases_n_iter=",
                     n_iter,
                     "_cutoff=",
                     cutoff,
                     "_FDR=",
                     FDR_logical,
                     "_permutation",
                     i,
                     ".rjob.sh"); cat(rjobsh, "\n")
    rjobConn <- file(paste0(dirname,"/",rjobsh))
    writeLines(c("#!/bin/sh -login",
                 "#SBATCH --mem=30GB",
                 paste0("#SBATCH --job-name=clust", i),
                 "#SBATCH --time=3:00:00",
                 "#SBATCH --nodes=1",
                 "#SBATCH --cpus-per-task=1",
                 "#SBATCH --account=wang-krishnan",
                 "",
                 "cd ../src",
                 "",
                 "module load GCC/8.3.0  OpenMPI/3.1.4  R/4.0.2",
                 "",
                 paste0("Rscript expandAndCluster_all.R ",
                        dg_path,
                        " ",
                        neighbors_path,
                        "/biogrid_permutation",
                        i,
                        "_neighbor_list.Rdata",
                        " ",
                        edgelist_path,
                        "/biogrid_permutation_",
                        i,
                        "_edgelist_1",
                        " ",
                        n_iter,
                        " ",
                        cutoff,
                        " ",
                        FDR_logical,
                        " ",
                        partition_type,
                        " ",
                        res,
                        " ",
                        results_path)),rjobConn)
    
    close(rjobConn)
    
    system(paste0("sbatch ", paste0(dirname,"/",rjobsh)))
    
    njobs <- system("squeue -u  mckimale | wc -l", intern=TRUE)
    njobs <- as.numeric(njobs)
    
    while(njobs > 1000) {
      Sys.sleep(360)
      njobs <- system("squeue -u mckimale | wc -l", intern=TRUE)
      njobs <- as.numeric(njobs)
    }
    print(getwd())
  }

