# This file will submit up to 1000 qsubs at once
library("tidyverse")

# path to folder with original disease gene list
#for args[1]
dg_path = "/mnt/research/compbio/krishnanlab/projects/chronic_inflammation/data/disease_gene_files/complexfaketraits"
setwd(dg_path)
dg_files = list.files(pattern = ".txt")

# path(s) to list of node neighbors from edgelist
# for args[2]
neighbors_path = "/mnt/research/compbio/krishnanlab/projects/chronic_inflammation/results/biogrid_noRWR/biogrid_neighbor_list.Rdata"
# path to edgelist
# for args[3]
edge_path = "/mnt/research/compbio/krishnanlab/projects/chronic_inflammation/data/biogrid/biogrid_entrez_edgelist.txt"
# results path for args[9]
results_path = "/mnt/research/compbio/krishnanlab/projects/chronic_inflammation/results/biogrid_noRWR"
# set pval cutoff
cutoff = 0.01
# use FDR?
FDR_logical = FALSE
# leiden alg partitian type
partition_type = "ModularityVertexPartition"
# resolution parameter
res = 0.1

setwd("/mnt/research/compbio/krishnanlab/projects/chronic_inflammation/run")

dirname <- "sbatches_expandAndCluster"
if(!dir.exists(dirname)){
  dir.create(dirname)
}

for(n_iter in c(1:5)) {
  for(file in dg_files){
    
    filename = sub(".txt", "", file)
    print(filename)
    rjobsh <- paste0(filename,
                     "_n_iter=",
                     n_iter,
                     "_cutoff=",
                     cutoff,
                     "_FDR=",
                     FDR_logical,
                     ".rjob.sh"); cat(rjobsh, "\n")
    rjobConn <- file(paste0(dirname,"/",rjobsh))
    writeLines(c("#!/bin/sh -login",
                 "#SBATCH --mem=100GB",
                 paste0("#SBATCH --job-name=", filename, n_iter),
                 "#SBATCH --time=0:30:00",
                 "#SBATCH --nodes=1",
                 "#SBATCH --cpus-per-task=1",
                 "#SBATCH --account=wang-krishnan",
                 "",
                 "cd /mnt/research/compbio/krishnanlab/projects/chronic_inflammation/src",
                 "",
                 "Rmodules",
                 "",
                 paste0("Rscript expandAndCluster.R ",
                        dg_path,
                        "/",
                        file,
                        " ",
                        neighbors_path,
                        " ",
                        edge_path,
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
    
    njobs <- system("squeue -u  hickeys6 | wc -l", intern=TRUE)
    njobs <- as.numeric(njobs)
    
    while(njobs > 1000) {
      Sys.sleep(360)
      njobs <- system("squeue -u hickeys6 | wc -l", intern=TRUE)
      njobs <- as.numeric(njobs)
    }
  }
}
