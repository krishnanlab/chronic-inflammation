############################################################
#' Creates 'seed gene files' from the uk biobank pascal results
#' @args[1] Uk biobank traits we are using
#' @args[2] directory with pascal output files
#' @args[3] Where disease_gene_files is located
############################################################
library(tidyverse)
args <- commandArgs(TRUE)

#convert to entrez
symbol2entrez = function(gene_symbols){
  require("org.Hs.eg.db")
  x <- org.Hs.egSYMBOL2EG
  included_symbols <- names(as.list(x))
  included_oi = gene_symbols[gene_symbols %in% included_symbols]
  xx <- as.list(x[included_oi])
  entrez_ids = unlist(xx)
  return(entrez_ids)
}

#Load file I made linking ids of interest and trait description
bio=read_tsv(as.character(args[1]))

#put underscores for every space in this description
#get rid of special characters of "-", "/" replacing them with _
#Hacky but it works for our files, and gets appropriate names
bio$X2=gsub("\\/","_",bio$X2)
bio$X2=gsub("\\-","_",bio$X2)
bio$X2=gsub(" ","_",bio$X2)
bio$X2=gsub(":","",bio$X2)
bio$X2=gsub("\\.","",bio$X2)
bio$X2=gsub("\\(","",bio$X2)
bio$X2=gsub("\\)","",bio$X2)
bio$X2=gsub("__","_",bio$X2)
bio$X2=gsub("__","_",bio$X2)

#Edit column 1 to have file name so I can get correct disease name
bio$X1=paste0(as.character(args[2]),"/",bio$X1,".gwas.imputed_v3.both_sexes.tsv.sum.genescores.txt")

pasfiles=list.files(as.character(args[2]),full.names=T)
#For each disease, get genes that have a pvalue < .001
for(file in pasfiles){
  pasresult=read_tsv(file,col_names=T)
  pasresult=pasresult %>% filter(pvalue<.001)
  pasentrez=symbol2entrez(pasresult$gene_symbol)
  dis=(filter(bio,X1==file))$X2
  pastrait=tibble(Gene=unname(pasentrez),Disease=dis)
  print(pastrait)
  write_tsv(pastrait,paste0(as.character(args[3]),"disease_gene_files/",dis,".txt"))
}
