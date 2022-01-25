# This file will submit up to 1000 qsubs at once
library("tidyverse")

# path to file with all random disease gene sets
#for args[1]
#random_path = "/mnt/research/compbio/krishnanlab/projects/chronic_inflammation/data/5000Expandedfaketraits_biogrid.tsv"
random_path = "../data_Zenodo/5000Expandedfaketraits_biogrid.tsv"
# load diseases of interest
# for args[2]
# object is sla1
load("../results/GenePlexus_parameters/diseases_with_SLA_models_cv_greater1.Rdata")

# path to graph
# for args[3]
edge_path = "../data_Zenodo/biogrid/BioGrid_igraph.Rdata"
netname = gsub("igraph.Rdata", "", basename(edge_path))

# results path for args[6]
results_path = "../results/GenePlexus_output/clusters_threshold.80"
# leiden alg partition type
partition_type = "ModularityVertexPartition"
# resolution parameter
res = 0.1

dirname <- "sbatches_clusterRandomGenes"
if(!dir.exists(dirname)){
  dir.create(dirname)
}

  for(disease in sla1[2:length(sla1)]){
    
    rjobsh <- paste0(disease,
                     ".rjob.sh"); cat(rjobsh, "\n")
    rjobConn <- file(paste0(dirname,"/",rjobsh))
    writeLines(c("#!/bin/sh -login",
                 "#SBATCH --mem=40GB",
                 paste0("#SBATCH --job-name=", disease, "_", netname),
                 paste0("#SBATCH --output=", dirname, "/", disease, "_", netname, ".out"),
                 "#SBATCH --time=1:00:00",
                 "#SBATCH --nodes=1",
                 "#SBATCH --cpus-per-task=1",
                 "",
                 "cd ../src",
                 "",
                 "ml -* GCC/8.3.0 OpenMPI/3.1.4 R/4.0.2",
                 "",
                 paste0("Rscript clusterRandomGenes.R ",
                        random_path,
                        " ",
                        disease,
                        " ",
                        edge_path,
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
  }
