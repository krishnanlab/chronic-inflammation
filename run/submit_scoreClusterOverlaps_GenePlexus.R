# This file will submit up to 1000 qsubs at once
library("tidyverse")
args <- commandArgs(TRUE)

dirname = "sbatches_ScoreClusterOverlaps_GenePlexus"
if(!dir.exists(dirname)){
  dir.create(dirname)}

clust_net = "ConsensusPathDB"
clust_net_dir = "ConsensusPathDB"
pred_net = "ConsensusPathDB"

# for args[5]
thresholds = 0.80

# for args[1]
cluster_dir = paste0("../results/prediction_clusters_same_graph/clusters/predicted_with", 
                     pred_net, 
                     "--clustered_on_", 
                     clust_net)

# for args[2]
ci_dir = "../results/GenePlexus_output"

ci_files = list.files(ci_dir,
                      full.names = TRUE,
                      pattern = paste0("--", pred_net ,"--Adjacency--GO--predictions.tsv"))

out_dir = "../results/prediction_clusters_same_graph"

# diseases of interest
# for args[4]

clust_files = list.files(cluster_dir)
diseases = gsub("--.*","",clust_files)
diseases = diseases[grep("Fake_", diseases, invert = TRUE)]

# args[6]
# genes in biogrid for args[6]
network_genes = paste0("../data/", clust_net_dir, "/", clust_net,"_genes.csv")

  for(ci in ci_files){
    
    cioi_name = sub(paste0("--", pred_net, "--Adjacency--GO--predictions.tsv"), "", basename(ci))
    
    for(thresh in thresholds){
      
      for(disease in diseases){
    
      score_file = paste0(out_dir, 
                          "/scores/",
                          cioi_name,
                          "_thresh=",
                          thresh,
                          "_predicted_with_",
                          pred_net,
                          "_clusteredOn_",
                          clust_net,
                          "/",
                          disease,
                          "_",
                          cioi_name,
                          "_thresh=",
                          thresh,
                          "_predicted_with_",
                          pred_net,
                          "_clusteredOn_",
                          clust_net,
                          "_enrichment_scores.csv")
      
      if(file.exists(score_file)){next}
  
  rjobsh <- paste0(cioi_name,
                   "_",
                   disease,
                   "_pred_",
                   pred_net,
                   "_clust_",
                   clust_net,
                   "_thresh=",
                   thresh,
                   ".rjob.sh"); cat(rjobsh, "\n")
  
  rjobConn <- file(paste0(dirname,"/",rjobsh))
  writeLines(c("#!/bin/sh -login",
               "#SBATCH --mem=100GB",
               paste0("#SBATCH --job-name=score_", disease),
               paste0("#SBATCH --output=", 
                      dirname, 
                      "/", 
                      disease, 
                      "_", 
                      cioi_name, 
                      "_pred_",
                      pred_net,
                      "_clust_",
                      clust_net, 
                      "_thresh=", 
                      thresh, 
                      ".out"),
               "#SBATCH --time=4:00:00",
               "#SBATCH --nodes=1",
               "#SBATCH --cpus-per-task=1",
               "#SBATCH --account=wang-krishnan",
               "",
               "cd ../src",
               "",
               "Rmodules",
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
  

  }
  }
}


