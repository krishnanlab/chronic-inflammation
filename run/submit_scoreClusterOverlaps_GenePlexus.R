# This file will submit up to 1000 qsubs at once
library("tidyverse")
args <- commandArgs(TRUE)

dirname = "sbatches_ScoreClusterOverlaps_GenePlexus"
if(!dir.exists(dirname)){
  dir.create(dirname)}

# for args[2]
ci = "../results/GenePlexus_output/predictions_cv_greater1/chronic_inflammation_go--STRING--Adjacency--GO--predictions.tsv"
# for args[5]
thresh = 0.80
# for args[1]
cluster_dir = "../results/GenePlexus_output/clusters_threshold.80"
# for args[3]
out_dir = "../results/GenePlexus_String_Adjacency"

# load diseases of interest
# for args[4]
# object is sla1
load("../results/GenePlexus_parameters/diseases_with_SLA_models_cv_greater1.Rdata")

# genes in biogrid for args[6]
network_genes = "../data_Zenodo/biogrid/BioGrid_genes.csv"

for(disease in sla1){
  
  cioi_name = sub("_clusters.csv", "", basename(ci))
  
  rjobsh <- paste0(cioi_name,
                   "_",
                   disease,
                   ".rjob.sh"); cat(rjobsh, "\n")
  
  rjobConn <- file(paste0(dirname,"/",rjobsh))
  writeLines(c("#!/bin/sh -login",
               "#SBATCH --mem=70GB",
               paste0("#SBATCH --job-name=score_", disease),
               paste0("#SBATCH --output=", dirname, "/", disease, "_", cioi_name, ".out"),
               "#SBATCH --time=3:59:00",
               "#SBATCH --nodes=1",
               "#SBATCH --cpus-per-task=1",
               "",
               "cd ../src",
               "",
               "ml -* GCC/8.3.0 OpenMPI/3.1.4 R/4.0.2",
               "",
               paste0("Rscript scoreClusterOverlaps_GenePlexus.R ",
                      cluster_dir,
                      " ",
                      ci,
                      " ",
                      out_dir,
                      " ",
                      disease,
                      " ",
                      thresh, 
                      " ",
                      network_genes)),rjobConn)
  
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


