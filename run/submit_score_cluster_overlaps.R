# This file will submit up to 1000 qsubs at once
#' @param args[1] path to chronic inflammation clusters for overlap
#' @param args[2] first job
#' @param args[3] last job

library("tidyverse")
args <- commandArgs(TRUE)

#for number of background genes
#change for w/e edgelist you're using
edge = read.delim("/mnt/research/compbio/krishnanlab/projects/chronic_inflammation/data/biogrid/biogrid_entrez_edgelist.txt", header = F)
bg = length(unique(c(edge$V1, edge$V2)))

dirname = "sbatches_ScoreClusterOverlaps"
if(!dir.exists(dirname)){
  dir.create(dirname)}

ci = "../results/biogrid_real_clusters/expanded_disease_genes_biogrid_iterations=2_FDR=FALSE_cutoff=0.01/clusters/biogrid_entrez_edgelist.txt_ModularityVertexPartition_res=NA/chronic_inflammation_go_clusters.csv"
# move real clusters to scratch with permuted clusters
# change cluster_dirs to scratch
cluster_dirs = dir("../../biogrid_hotnet2_permutations/biogrid_hotnet2_permutations_noRWR/expanded_disease_genes_biogrid_iterations=2_FDR=FALSE_cutoff=0.01/clusters",
                   full.names = TRUE)
out_dir = "../results/overlaps"

#first = as.numeric(args[2])
#last = as.numeric(args[3])

#cluster_dirs = cluster_dirs[first:last]
cluster_dirs=cluster_dirs[1:1000]
  count = 0
  for(clust in cluster_dirs){
      count = count + 1
      cioi_name = sub("_clusters.csv", "", basename(ci))
      clust_dir_name = sub("_ModularityVertexPartition_res=NA", "", basename(clust))
     
      rjobsh <- paste0(cioi_name,
                       "_",
                       clust_dir_name,
                       ".rjob.sh"); cat(rjobsh, "\n")
      
      rjobConn <- file(paste0(dirname,"/",rjobsh))
      writeLines(c("#!/bin/sh -login",
                   "#SBATCH --mem=30GB",
                   paste0("#SBATCH --job-name=score", count),
                   "#SBATCH --time=1:00:00",
                   "#SBATCH --nodes=1",
                   "#SBATCH --cpus-per-task=1",
                   "#SBATCH --account=wang-krishnan",
                   "",
                   "cd ../src",
                   "",
                   "module load GCC/8.3.0  OpenMPI/3.1.4  R/4.0.2",
                   "",
                   paste0("Rscript score_cluster_overlaps.R ",
                          clust,
                          " ",
                          ci,
                          " ",
                          out_dir,
                          " ",
                          bg)),rjobConn)
      
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


