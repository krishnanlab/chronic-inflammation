# This file will submit up to 1000 qsubs at once
library("tidyverse")
# number of hypergeo iterations
setwd("../run")

dirname <- "sbatches_GOBP_enriched_clusters_GenePlexus"
if(!dir.exists(dirname)){
  dir.create(dirname)
}

netname = "STRING-EXP"
netdir = "string-exp"
pred_net = "STRING-EXP"

# path to background genes
background = paste0("../data/", 
                    netdir, 
                    "/", 
                    netname,
                    "_genes.csv")

#outdir
outdir = "../results/prediction_clusters_same_graph/GOBP_enrichment"

cluster_dir = paste0("../results/prediction_clusters_same_graph/clusters/predicted_with", 
                     pred_net,
                     "--clustered_on_",
                     netname)

cluster_files = list.files(cluster_dir, full.names = TRUE)
cluster_files = cluster_files[grep("Fake_", cluster_files, invert=TRUE)]
redo = c("Own_or_rent_accommodation_lived_in_Rent_from_local_authority,_local_council,_housing_association", "Gas_or_solid_fuel_cooking_heating_An_open_solid_fuel_fire_that_you_use_regularly_in_winter_time", "Malignant_neoplasm_of_pancreas", "PARKINSON_DISEASE_2_AUTOSOMAL_RECESSIVE_JUVENILE")

cluster_files = cluster_files[grep(paste(redo,collapse="|"), cluster_files)]


  for(file in cluster_files){
    
    filename = sub(".csv", "", basename(file))
    print(filename)
    rjobsh <- paste0(filename,
                     "_GOBP.rjob.sh"); cat(rjobsh, "\n")
    rjobConn <- file(paste0(dirname,"/",rjobsh))
    writeLines(c("#!/bin/sh -login",
                 "#SBATCH --mem=128GB",
                 paste0("#SBATCH --job-name=", filename),
                 paste0("#SBATCH --output=", dirname, "/", filename, "_GOBP.out"),
                 "#SBATCH --time=3:00:00",
                 "#SBATCH --nodes=1",
                 "#SBATCH --cpus-per-task=1",
                 "#SBATCH --account=wang-krishnan",
                 "",
                 "cd ../src",
                 "",
                 "Rmodules",
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

  }

