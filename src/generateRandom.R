#for each trait, generate 1000 fake traits usng same node degree
library(tidyverse)
library(parallel)
library(igraph)
library(org.Hs.eg.db)
source("chronic_inflammation_functions.R")

generatetraits<-function(disease,origseeds){
  set.seed(35)
  seeds=origseeds %>% filter(Disease==disease)
  seeds$Entrez=as.character(seeds$Gene)
  seeds=seeds %>%dplyr::select(-Gene)
  print(seeds)
  #strings=read_tsv("/mnt/research/compbio/krishnanlab/projects/chronic_inflammation/data/string/biogrid_gene_index.txt",col_names=F) %>% dplyr::select(-X1)
  #string = read_tsv("/mnt/research/compbio/krishnanlab/data/ukb-gwas_neallab/STRING.edg",col_types="ccd",col_names=F)
  string = read_tsv("../data_Zenodo/biogrid/biogrid_entrez_edgelist.txt",col_types="cc",col_names=F)
  
  stringgenes= c(string$X1,string$X2) %>% unique %>% as.character %>% as_tibble
  colnames(stringgenes)="Entrez"
  #colnames(strings)="Entrez"
  #strings$Entrez=as.character(strings$Entrez)
  #Get genes that are not seed genes but in string
  syms=entrez2symbol(stringgenes$Entrez)
  allstringentrez=names(syms)
  
  #genes in string not in seeds
  stringnonseeds=allstringentrez[!(allstringentrez %in% seeds$Entrez)]
  
  
  #node degree is how many connections it has
  edgelist=string %>% dplyr::select(X1,X2)
  g=graph_from_edgelist(edgelist %>% as.matrix)
  degreesnamedvector=degree(g,V(g))
  degreestibble=as.data.frame(degreesnamedvector) 
  degreestibble$Entrez=rownames(degreestibble)
  degreestibble= degreestibble %>% as_tibble %>% dplyr::select(Entrez,degreesnamedvector)
  
  degreestodraw=left_join(seeds,degreestibble,by="Entrez")
  colnames(degreestodraw)[3]="degreeneeded"
  #Gets rid of seeds from drawng pool
  #drawgenesfrom=filter(degreestibble,!(Entrez%in%degreestodraw$Entrez))
  #keep seeds
  drawgenesfrom=degreestibble
  answer=tibble()
  loop<-function(i){
    faketraitgenes=character()
    #for each gene, find a replacement gene
    for(j in 1:nrow(degreestodraw)){
      #current degree to find a gene for
      deg=pull(degreestodraw[j,"degreeneeded"])
      #seed gene not in string
      if(is.na(deg)){
        next
      }
      #genes that have this degree that I can draw from
      #Get degreesnamedvector that are within a range if <3 options
      
      ens=filter(drawgenesfrom,degreesnamedvector==deg)$Entrez
      #print(length(ens))
      while(length(ens)<3){
        #print(j)
        #print(length(ens))
        #print(ens)
        #sort to no longer deal with indices, aka whatever is next to the min
        #will be genes I am interested in
        #If want to get closest degrees, then put %>% unique after this line
        #otherwise we get the closest nearby genes
        de=sort(drawgenesfrom$degreesnamedvector)
        closestindex=which.min(abs(de - deg))
        #get 5 closest node degree indices
        maxnum=closestindex+30
        if(maxnum>length(de)){
          maxnum=length(de)
        }
        minnum=closestindex-30
        if(minnum<1){
          minnum=1
        }
        degoptions=de[minnum:maxnum]
        
        ens=filter(drawgenesfrom,degreesnamedvector %in% unique(degoptions))$Entrez
        #print(length(ens))
      }
      #If no gene, need to get closest one
      while(length(ens)==0){
        closestdegindex=which.min(abs(drawgenesfrom$degreesnamedvector - deg))
        closestdeg=drawgenesfrom$degreesnamedvector[closestdegindex]
        ens=filter(drawgenesfrom,degreesnamedvector==closestdeg)$Entrez
      }
      en=sample(ens,1)
      faketraitgenes=c(faketraitgenes,en)
      
    }
    #print(nrow(genes))
    #distinct because if duplicate because closest gene is to different high degree seed genes
    iter=tibble(Gene=faketraitgenes,Disease=paste0("Fake_",disease,"_",i)) %>% distinct
    #print(tibble(Gene=faketraitgenes,Disease=paste0("Fake_",disease,"_",i)))
    #print(i)
   return(iter) 
  }
  is=1:5000
  iters=mclapply(is,loop,mc.cores=detectCores()-1)
  answer=bind_rows(iters)
  print(answer)
  return(answer)
}
fi=list.files("../data/disease_gene_files",full.names=T)
origseeds=mclapply(fi,read_tsv) %>% bind_rows


diseases=unique(origseeds$Disease)
answer=lapply(diseases,generatetraits,origseeds)
answers=bind_rows(answer)
write_tsv(answers,"../data/5000Expandedfaketraits_biogrid.tsv")