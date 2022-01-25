# This file will submit up to 1000 qsubs at once
library("tidyverse")
# number of hypergeo iterations
dirname <- "sbatches_GOBP_enriched_clusters_GenePlexus"
if(!dir.exists(dirname)){
  dir.create(dirname)
}

# path to background genes
#background = "/mnt/research/compbio/krishnanlab/projects/chronic_inflammation/data/biogrid/BioGrid_genes.csv"
background = "../data_Zenodo/biogrid/BioGrid_genes.csv"

#outdir
outdir = "../results/GenePlexus_output/GOBP_enrichment"

cluster_dir = "../results/GenePlexus_output/clusters_threshold.80"
cluster_files = list.files(cluster_dir, full.names = TRUE, pattern = ".csv")
  
  for(file in cluster_files){
    
    filename = sub(".csv", "", basename(file))
    print(filename)
    rjobsh <- paste0(filename,
                     "GOBP.rjob.sh"); cat(rjobsh, "\n")
    rjobConn <- file(paste0(dirname,"/",rjobsh))
    writeLines(c("#!/bin/sh -login",
                 "#SBATCH --mem=60GB",
                 paste0("#SBATCH --job-name=", filename),
                 paste0("#SBATCH --output=", dirname, "/", filename, "_GOBP.out"),
                 "#SBATCH --time=2:00:00",
                 "#SBATCH --nodes=1",
                 "#SBATCH --cpus-per-task=1",
                 "",
                 "cd ../src",
                 "",
                 "ml -* GCC/8.3.0 OpenMPI/3.1.4 R/4.0.2",
                 "",
                 paste0("Rscript find_GOBP_enriched_clusters_GenePlexus.R ",
                        file,
                        " ",
                        background,
                        " ",
                        outdir
                 )),rjobConn)
    
    close(rjobConn)
    
    system(paste0("sbatch ", paste0(dirname,"/",rjobsh)))
    
    njobs <- system("squeue -u mckimale | wc -l", intern=TRUE)
    njobs <- as.numeric(njobs)
    
    while(njobs > 1000) {
      Sys.sleep(360)
      njobs <- system("squeue -u mckimale | wc -l", intern=TRUE)
      njobs <- as.numeric(njobs)
    }
  }

