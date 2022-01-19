# This file will submit up to 1000 qsubs at once
library("tidyverse")

# path to folder with gene plexus predictions
#for args[1]
gp_path = "../results/GenePlexus_output/predictions_cv_greater1"
#setwd(gp_path)
#gp_files = list.files(pattern = ".tsv")
gp_files = list.files(gp_path,pattern=".tsv")

# path(s) to list of node neighbors from edgelist
# for args[2]
#thresholds = c("mccf1", seq(0.5, 0.95, .05))
#Using .8 because thats what we ended up using, origianlly ran with sequence
thresholds=c("0.8",)

# path to graph
# for args[3]
edge_path = "../Zenodo/biogrid/BioGrid_igraph.Rdata"
netname = gsub("igraph.Rdata", "", basename(edge_path))
# results path for args[6]
results_path ="../results/GenePlexus_output/"
# leiden alg partition type
partition_type = "ModularityVertexPartition"
# resolution parameter
res = 0.1

#setwd("../run")
print(getwd())
dirname <- "sbatches_filterAndClusterGeneplexusPredictions"
if(!dir.exists(dirname)){
  dir.create(dirname)
}

for(thresh in thresholds) {
  for(file in gp_files){
    
    filename = sub("--predictions.tsv", "", file)
    rjobsh <- paste0(filename,
                     "_treshold=",
                     thresh,
                     ".rjob.sh"); cat(rjobsh, "\n")
    rjobConn <- file(paste0(dirname,"/",rjobsh))
    writeLines(c("#!/bin/sh -login",
                 "#SBATCH --mem=40GB",
                 paste0("#SBATCH --job-name=", filename, "_", thresh, "_", netname),
                 paste0("#SBATCH --output=", dirname, "/", filename, "_", thresh, "_", netname, ".out"),
                 "#SBATCH --time=0:30:00",
                 "#SBATCH --nodes=1",
                 "#SBATCH --cpus-per-task=1",
                 "#SBATCH --exclude=qml-000,qml-002,qml-003",
                 "#SBATCH --account=wang-krishnan",
                 "",
                 "cd ../src",
                 "",
                 "ml -* GCC/8.3.0 OpenMPI/3.1.4 R/4.0.2",
                 "",
                 paste0("Rscript filterAndClusterGeneplexusPredications.R ",
                        gp_path,
                        "/",
                        file,
                        " ",
                        thresh,
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
}
