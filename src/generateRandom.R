#for each trait, generate 1000 fake traits usng same node degree
library(tidyverse)
library(parallel)
library(igraph)
library(org.Hs.eg.db)
source("chronic_inflammation_functions.R")

generatetraits<-function(disease,origseeds,degreestibble,toDraw){
  set.seed(35)
  seeds=origseeds %>% filter(Disease==disease)
  seeds$Entrez=as.character(seeds$Gene)
  seeds=seeds %>%dplyr::select(-Gene)
  print(seeds)

  #keep only our seed genes
  drawset=filter(toDraw,ToReplace %in% seeds$Entrez)
  
  #degreestodraw=left_join(seeds,degreestibble,by="Entrez")
  #colnames(degreestodraw)[3]="degreeneeded"
  
  answer=tibble()
  loop<-function(i){
    #print(i)
    sampled=drawset %>% group_by(ToReplace) %>% slice_sample(n=1)
    #print(nrow(genes))
    #distinct because sometimes I will get the same gene multiple times in drawing
    iter=tibble(Gene=sampled$BinnedGenes,Disease=paste0("Fake_",disease,"_",i)) %>% distinct
    iter
    #print(tibble(Gene=faketraitgenes,Disease=paste0("Fake_",disease,"_",i)))
    #print(i)
   return(iter) 
  }
  is=1:5000
  iters=mclapply(is,loop,mc.cores=detectCores())
  answer=bind_rows(iters)
  print(answer)
  return(answer)
}
#seeds=read_csv("/mnt/research/compbio/krishnanlab/projects/chronic_inflammation/data/all_seed_genes.csv")
#fi=list.files("/mnt/research/compbio/krishnanlab/projects/chronic_inflammation/data/disease_gene_files",full.names=T)
fi=list.files("/mnt/research/compbio/krishnanlab/projects/chronic_inflammation/data/disease_gene_files_SLA_thresh_0.80",full.names=T)
fi=fi[fi!="/mnt/research/compbio/krishnanlab/projects/chronic_inflammation/data/disease_gene_files/Random_biogrid_matched"]
origseeds=mclapply(fi,read_tsv) %>% bind_rows


diseases=unique(origseeds$Disease)
#diseases=diseases[!(diseases %in% c("Banana_Intake","chronic_inflammation_gene_shot_pubs_greater10","chronic_inflammation_gene_shot",
#                                    "chronic_inflammation_go","Driving_Speed","East_UK","Equal_Hand_Phone","Faketrait_1","Faketrait_2","Faketrait_3","Faketrait_4","Faketrait_5",
#                                    "Faketrait_6","Faketrait_7","Faketrait_8","Faketrait_9","Height","Left_Hand_Phone","North_UK","Right_Hand_Phone"))]
#degreestibble=tibble()
#getdegrees=function(gene){
#  fi=filter(string,X1==gene | X2==gene)
#  degree=nrow(fi)
#  avgweight=mean(fi$X3)
#  degreestibble=tibble(Entrez=gene,Degree=degree,Average_Weight=avgweight)
#  print(gene)
#  return(degreestibble)
#}

#Choose network
#string = read_tsv("../data/string/STRING_human_entrez_edgelist.txt",col_types="ccd",col_names=F)
#string = read_tsv("../data/string-exp/STRING-EXP_human_entrez_edgelist.txt",col_types="ccd",col_names=F)
#string = read_tsv("../data/ConsensusPathDB/ConsensusPathDB_human_entrez_edgelist.txt",col_types="ccd",col_names=F)
string = read_tsv("../data/biogrid/biogrid_entrez_edgelist.txt",col_types="ccd",col_names=F)

stringgenes= c(string$X1,string$X2) %>% unique %>% as.character %>% as_tibble
colnames(stringgenes)="Entrez"
#colnames(strings)="Entrez"
#strings$Entrez=as.character(strings$Entrez)
#Get genes that are not seed genes but in string
syms=entrez2symbol(stringgenes$Entrez)
allstringentrez=names(syms)

#node degree is how many connections it has
edgelist=string %>% dplyr::select(X1,X2)
g=graph_from_edgelist(edgelist %>% as.matrix)
degreesnamedvector=degree(g,V(g))
degreestibble=as.data.frame(degreesnamedvector) 
degreestibble$Entrez=rownames(degreestibble)
degreestibble= degreestibble %>% as_tibble %>% dplyr::select(Entrez,degreesnamedvector)

#node degree, and how many times it appears
#Make a tibble where I bin up each gene that
#has a low node degree before this loop. Then give that gene a list I can
#draw from. This will speed things up tremendously
tab=table(degreestibble$degreesnamedvector)
degappear=as.data.frame(tab) %>% as_tibble
colnames(degappear)=c("Degree","Freq")
degappear$Degree=as.numeric(as.character(degappear$Degree))

colnames(degreestibble)[2]="Degree"
degreestibble = degreestibble %>% left_join(degappear,by="Degree")
#make a tibble where it will be:
#Gene I want to replace | Possible gene to replace it with
#Where in the loop I call mclapply for, I filter this tibble for the gene
#I want to replace, and then draw randomly from column 2


makeBins<-function(i){
  #print(i)
  #If the number of genes with this degree is less than, then we get more genes
  fin=tibble()
  hold=degreestibble
  degreestibble= degreestibble %>% arrange(Degree)
  curgene=pull(degreestibble[i,1])
  deg=pull(degreestibble[i,2])
  freq=pull(degreestibble[i,3])
  if(freq<10){
    #If I hae very few, then add more. If only 1 gene, then add more options
    #So every gene should have at least 10 genes in the pool to sample from
    #delete this line and restore the commented out lines (and get rid of the lines
    #that replace them that have the toadd var in them) to get old script back
    toadd=10/2
    
    de=sort(degreestibble$Degree)
    closestindex=which.min(abs(de - deg))

    #maxnum=closestindex+ (2*freq)
    maxnum=closestindex + toadd -1 #Add a -1 to get only 10 rather than 11
    addtomin=0
    if(maxnum>length(de)){
      maxnum=length(de)
      addtomin=toadd
      #addtomin=2*freq
    }
    #minnum=closestindex- (2*freq) - addtomin
    minnum=closestindex- toadd - addtomin
    if(minnum<1){
      minnum=1
      #maxnum=maxnum+ (2*freq)
      maxnum=maxnum+ (toadd)
    }
    #genes to get
    curselect=degreestibble[minnum:maxnum,]
    fin=bind_rows(fin,tibble(ToReplace=rep(curgene,nrow(curselect)),BinnedGenes=curselect$Entrez))
  }else{
      #Just add this curgene, and the other genes with taht same degree as binned
      subs=filter(degreestibble,Degree==deg)
      fin=bind_rows(fin,tibble(ToReplace=rep(curgene,nrow(subs)),BinnedGenes=subs$Entrez))
      
  }
  return(fin)
}
i=1:nrow(degreestibble)
an=mclapply(1:nrow(degreestibble),makeBins,mc.cores=detectCores())
fin=bind_rows(an)

#These 2 lines will let you see if I got the genes I needed or not. This is an example
#where enough genes so should just make a bin of genes with same degree
fin %>% filter(ToReplace==1)
filter(degreestibble,Degree==396)

toDraw=fin

answer=lapply(diseases,generatetraits,origseeds,degreestibble=degreestibble,toDraw=fin)
answers=bind_rows(answer)

#write_tsv(answers,"../data/5000Expandedfaketraits_STRING.tsv")
#write_tsv(answers,"../data/5000Expandedfaketraits_STRING-EXP.tsv")
#write_tsv(answers,"../data/5000Expandedfaketraits_ConsensusPathDB.tsv")
write_tsv(answers,"../data/5000Expandedfaketraits_biogrid.tsv")