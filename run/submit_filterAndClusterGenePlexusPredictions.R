# This file will submit up to 1000 qsubs at once
library("tidyverse")

# geneplexus prediction network
pred_net = "STRING"
clust_dir = "string"
clust_net = "STRING"
# weighted network? for args[7]
weight = T

# path to folder with gene plexus predictions
#for args[1]
gp_path = "/mnt/research/compbio/krishnanlab/projects/chronic_inflammation/results/GenePlexus_output"
setwd(gp_path)
gp_files = list.files(pattern = ".tsv")
gp_files = gp_files[grep(paste0(pred_net, "--"), gp_files)]
gp_files = gp_files[grep("inflamm", gp_files, invert = T)]

# path(s) to list of node neighbors from edgelist
# for args[2]
# thresholds = c("mccf1", seq(0.5, 0.95, .05))
thresholds = .80
# path to graph
# for args[3]
edge_path = paste0("/mnt/research/compbio/krishnanlab/projects/chronic_inflammation/data/",
                   clust_dir, 
                   "/",
                   clust_net,
                   "_igraph.Rdata")

# results path for args[6]
results_path = "/mnt/research/compbio/krishnanlab/projects/chronic_inflammation/results/prediction_clusters_same_graph"
# leiden alg partition type
partition_type = "ModularityVertexPartition"
# resolution parameter
res = 0.1

setwd("/mnt/research/compbio/krishnanlab/projects/chronic_inflammation/run")

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
                     "_cluster_on_",
                     clust_net,
                     ".rjob.sh"); cat(rjobsh, "\n")
    rjobConn <- file(paste0(dirname,"/",rjobsh))
    writeLines(c("#!/bin/sh -login",
                 "#SBATCH --mem=100GB",
                 paste0("#SBATCH --job-name=", filename, "_", thresh, "_", clust_net),
                 paste0("#SBATCH --output=", dirname, "/", filename, "_", thresh, "_", clust_net, ".out"),
                 "#SBATCH --time=3:00:00",
                 "#SBATCH --nodes=1",
                 "#SBATCH --cpus-per-task=1",
                 "#SBATCH --account=wang-krishnan",
                 "",
                 "cd /mnt/research/compbio/krishnanlab/projects/chronic_inflammation/src",
                 "",
                 "Rmodules",
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
                        results_path,
                        " ",
                        weight)),rjobConn)
    
    close(rjobConn)
    
    system(paste0("sbatch ", paste0(dirname,"/",rjobsh)))

  }
}
