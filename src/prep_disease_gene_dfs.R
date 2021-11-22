############################################################
#' For each disease ID outputs a tab delimted file 
#' with two columns (Gene and Disease) to a folder
#' named "disease_gene_files"
#' @args[1] a text file with one column (no header)
#' containing the disease ids of interest from disgenet
#' @args[2] dir where "disease_gene_files" folder will go
############################################################

library(tidyverse)
args <- commandArgs(TRUE)
source("chronic_inflammation_functions.R")
# load disgenet table
dg = read.delim("../data/disgenet_2020-10_curated_gene_disease_associations.tsv")
# load disease ids of interest
diseases = read.delim(args[1], header = F)
diseases = diseases$V1
# set up the output_dir
output_dir = args[2]
if (!dir.exists(output_dir)) {dir.create(output_dir)}

disease_dir = paste0(output_dir, "/disease_gene_files")
if (!dir.exists(disease_dir)) {dir.create(disease_dir)}

for(d in diseases){
  
  filt = dg %>%
    filter(diseaseId == d) 
  
  disease_name = unique(filt$diseaseName)
  disease_name = gsub(" ", "_", disease_name)
  disease_name = gsub("([_])|[[:punct:]]", "\\1", disease_name)
  
  df = filt %>%
    select(geneId)
  
  colnames(df) = "Gene"
  df$Disease = disease_name
  
  write.table(df, 
              file = paste0(disease_dir, "/", disease_name, ".txt"),
              sep = "\t",
              row.names = F,
              quote = F)
}

# chronic inflammation seed genes
# chronic_inflammation_geneshot_npubs_greater10
cigs10 = read.delim("../data/chronic_inflammation_genes/chronic_inflammation_geneshot_npubs_greater10.tsv")
gene = cigs10$Gene
entrez = symbol2entrez(gene)
cigs10 = data.frame(Gene = entrez,
                    Disease = "chronic_inflammation_gene_shot_pubs_greater10")
write.table(cigs10, 
            file = paste0(disease_dir, "/", "chronic_inflammation_gene_shot_pubs_greater10.txt"),
            sep = "\t",
            row.names = F,
            quote = F)
# chronic-inflammation_geneshot-generif_genes
cigs = read.delim("../data/chronic_inflammation_genes/chronic-inflammation_geneshot-generif_genes.tsv")
gene = cigs$Gene
entrez = symbol2entrez(gene)
cigs = data.frame(Gene = entrez,
                    Disease = "chronic_inflammation_gene_shot")
write.table(cigs, 
            file = paste0(disease_dir, "/", "chronic_inflammation_gene_shot.txt"),
            sep = "\t",
            row.names = F,
            quote = F)

# human_chronic_inflammatory_reponse_GO
cigo = read.delim("../data/chronic_inflammation_genes/human_chronic_inflammatory_reponse_GO.txt",header=F)
gene = cigo$V2
entrez = symbol2entrez(gene)
cigo = data.frame(Gene = entrez,
                  Disease = "chronic_inflammation_go")
write.table(cigo, 
            file = paste0(disease_dir, "/", "chronic_inflammation_go.txt"),
            sep = "\t",
            row.names = F,
            quote = F)

