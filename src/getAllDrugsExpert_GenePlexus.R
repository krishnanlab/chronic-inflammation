#Will e the new file I use to get disease/inflam drugs, as well as using synonyms to search dgidb for interactions

#This function take each drugcentral drug, get synonyms, and see if the synonym is in
#expertcurated. If it is, overwrite it in expertcurated with the preferred name I am
#using in drugcentral
getsynonyms<-function(curid){
  #the id is a struct id from drugcentral, for that id, get all instances in
  #synonyms.csv
  print(curid)
  namesfordrug=(synonyms %>% filter(id==curid))$name
  foundnames=expertcurated %>% filter(drugName %in% namesfordrug)
  #take foundnames, replace each row with name I use in drucentral, and then
  #return this
  foundnames$drugName=unique((drugcentral %>% filter(STRUCT_ID==curid))$DRUG_NAME)
  return(foundnames)
}
#Get for each row the num of cols greater than
getpvalues<-function(row,faketraitnum){
  
  #number of permutations the real score is less than
  #get less than so I can just divide. If less than all, I will get 1000. 1000/1000=1
  #I want cols that are less than score, because that means the score I have is better
  #, so >= beucase if equal, then it will get more cols and lower the pval
  cols_greaterequal_score=pivoted[row,2:ncol(pivoted)] %>% select_if(~any(.>=pull(pivoted[row,2])))
    
  #divide by 1000 here. Why? Because if I lose permutations because no clusters,
  #I should count them since the only reason I lack them is because no clusters were found,
  #thus no drugs were found
  #add 1 to 1000 for the real
  pvalue=as.double(ncol(cols_greaterequal_score) / (faketraitnum +1 ))
  #pvalue=as.double(ncol(cols_score_less_than) / (ncol(pivoted) - 2))
  #print(pvalue)
  return(pvalue)
}

#Take in clustermatrix, and for each permutation I get the drugs that hit each disease
#and the score. 
#druglist because ill just have every drug possible, and after this delete rows
#where all NA for a disease
getdiseasedrugs<-function(graphnum){
  whattoget=unique(diseaseclusters$Graph)[graphnum]
  
  #if real, only get clusters from diseaseclusters that are in the names of the clustermatrix
  #else get all clusters
  if(whattoget=="biogrid"){
    curgraph=diseaseclusters %>% filter(Graph==whattoget,Cluster %in% clustermatrix$Disease_Cluster)
  }else{
    curgraph=diseaseclusters %>% filter(Graph==whattoget)  
  }
  #set up tibble
  toreturn=tibble()
  #for each disease in curgraph
  for(disease in unique(curgraph$Disease)){
    #print(disease)
    #filter for the disease
    curdis=filter(curgraph,Disease==disease)
    #get all drugs that are in the targets for this disease
    founddrugs=drugcentral %>% filter(Entrez %in% curdis$Gene)
    colnames(curdis)[1]="Entrez"
    #Get cluster by joining found drugs with curdis
    founddrugs=founddrugs %>% left_join(curdis,by="Entrez")
    
    #For the drugs that are found, find targets that appear. Don't need to filter
    #again, just see how many rows the drug is in drugcenral, and how many rows
    #it is in founddrugs
    drugs=character()
    scores=numeric()
    entrezfound=character()
    symbolfound=character()
    Clusterforgenes=character()
    for(drug in unique(founddrugs$DRUG_NAME)){
      drugs=c(drugs,drug)
      scores=c(scores,length(founddrugs$DRUG_NAME[founddrugs$DRUG_NAME==drug]) / length(drugcentral$DRUG_NAME[drugcentral$DRUG_NAME==drug]))
      entrezfound=c(entrezfound,paste0(founddrugs$Entrez[founddrugs$DRUG_NAME==drug],collapse=","))
      symbolfound=c(symbolfound,paste0(founddrugs$GENE[founddrugs$DRUG_NAME==drug],collapse=","))
      Clusterforgenes=c(Clusterforgenes,paste0(founddrugs$Cluster[founddrugs$DRUG_NAME==drug],collapse=","))
    }
    toreturn=rbind(toreturn,tibble(Network=whattoget,Disease=disease,Drugs=drugs,Scores=scores,EntrezFound=entrezfound,SymbolFound=symbolfound,
                                   Clustersforgenes=Clusterforgenes))
    
  }
  print(toreturn)
  
  #Return large matrix where disease |  drug | col of permuted/real | score
  #I will cbind all of the results in the list
  return(toreturn)
}

#Get drugs for all traits, real or fake 
getDiseaseFromFake<-function(trait){
  #subset for fake trait
  curtrait=diseaseclusters %>% filter(Disease==trait)
  #if a real trait, further subset based on clustermatrix$Disease_Cluster
  if(grepl("Fake",trait)==F){
    print(trait)
    print(nrow(curtrait))
    curtrait = curtrait %>% filter(Cluster %in% clustermatrix$Disease_Cluster)
    print(nrow(curtrait))
  }
  
  #search for drugs from targets of fake trait
  founddrugs = drugcentral %>% filter(Entrez %in% curtrait$Gene)
  colnames(curtrait)[1]="Entrez"
  founddrugs = left_join(founddrugs,curtrait,by="Entrez")
  toreturn=tibble()
    #curdis=diseaseclusters %>% filter(Disease==i)
    
    drugs=character()
    scores=numeric()
    entrezfound=character()
    symbolfound=character()
    Clusterforgenes=character()
    for(drug in unique(founddrugs$DRUG_NAME)){
      drugs=c(drugs,drug)
      scores=c(scores,length(founddrugs$DRUG_NAME[founddrugs$DRUG_NAME==drug]) / length(drugcentral$DRUG_NAME[drugcentral$DRUG_NAME==drug]))
      entrezfound=c(entrezfound,paste0(founddrugs$Entrez[founddrugs$DRUG_NAME==drug],collapse="|"))
      symbolfound=c(symbolfound,paste0(founddrugs$GENE[founddrugs$DRUG_NAME==drug],collapse="|"))
      Clusterforgenes=c(Clusterforgenes,paste0(founddrugs$Cluster[founddrugs$DRUG_NAME==drug],collapse=","))
    }
    toreturn=rbind(toreturn,tibble(Trait=trait,Drugs=drugs,Scores=scores,EntrezFound=entrezfound,SymbolFound=symbolfound,
                                   Clustersforgenes=Clusterforgenes))
  #print(trait)
  return(toreturn)
}



library(tidyverse)
library(parallel)
#inflamclusters=read_csv("/mnt/research/compbio/krishnanlab/projects/chronic_inflammation/results/biogrid_noRWR/expanded_disease_genes_biogrid_iterations=2_FDR=FALSE_cutoff=0.01/clusters/chronic_inflammation_go_clusters.csv")%>%
#  select(Cluster,Gene)
scores=read_tsv("../output/diseasescorevector.tsv",col_names=T)
#load expertcurated
expertcurated=read_tsv("./../data/dgidb/expertcurated_entrez.tsv")
expertcurated = expertcurated %>% dplyr::select(drugName,Entrez,expert) %>% distinct()
drugcentral=read_tsv("../data/drugcentral/drugcentral_entrez.tsv")

synonyms=read_csv("../data/drugcentral/databasecsvs/synonyms.csv")
synonyms$name=toupper(synonyms$name)
synonyms=filter(synonyms,id %in% drugcentral$STRUCT_ID)
#number of ids from drugcentral
ids=unique(drugcentral$STRUCT_ID)
out=mclapply(ids,getsynonyms,mc.cores=8)
expertcurated=bind_rows(out) %>% distinct
colnames(expertcurated)=c("DRUG_NAME","Entrez","expert")
#filter out nonexpert
drugcentral=drugcentral %>% left_join(expertcurated,by=c("DRUG_NAME","Entrez")) %>% filter(!is.na(expert)) %>% select(-expert,-TDL) %>% distinct
drugcentral$Entrez=as.character(drugcentral$Entrez)
#need to find drugs for each real and permuted networks
#for each disease, get drug, see how many targets lie in the significant'y overlapping
#disease clusters, and get score based on percent of targets that are in the clusters for that
#disease

#load("/mnt/research/compbio/krishnanlab/projects/chronic_inflammation/results/all_cluster_assignments_for_alex.Rdata")
#Do this stuff if Stephanie did not filter clusters out
all_clusters_df=read_csv("/mnt/research/compbio/krishnanlab/projects/chronic_inflammation/results/GenePlexus_String_Adjacency/relevant_gene_cluster_assigments.csv")
diseaseclusters=all_clusters_df %>% as_tibble %>% filter(Disease!="chronic_inflammation_go") %>% select(-X1,-ClusterGraph)
diseaseclusters$Gene=as.character(diseaseclusters$Gene)
#rm(all_clusters_df)
#tokeep=diseaseclusters %>% count(Graph,Cluster) %>% filter(n>=5) %>% select(-n)
#diseaseclusters=tokeep %>% left_join(diseaseclusters,by=c("Graph","Cluster")) %>% select(Gene,Disease,Cluster,ClusterID,Method,Graph,Iterations,expandFDR,expandCutoff,Resolution,Control)
#rm(tokeep)


clustermatrix=read_csv("/mnt/research/compbio/krishnanlab/projects/chronic_inflammation/results/GenePlexus_String_Adjacency/final_for_alex.csv")#%>% filter(phyperFDR<.05) %>% filter(PermutedFDR<.05)
faketraits=diseaseclusters %>% filter(grepl("Fake",Disease))
#adding real clusters to above, also filter for only relevant clusters from clustermatrix
realtraits=diseaseclusters %>% filter(!grepl("Fake",Disease)) %>% filter(Cluster %in%clustermatrix$Disease_Cluster)
toget=c(unique(realtraits$Disease),unique(faketraits$Disease))

out=mclapply(toget,getDiseaseFromFake,mc.cores=14)
#each thing in out is each permutation, so I want to rbind everything
#1-27 are real traits
drugscores=bind_rows(out)
realtraitscores=filter(drugscores,!grepl("Fake",Trait))
faketraitscores=filter(drugscores,grepl("Fake",Trait))
#Loop for dataframe where I get each disease that is in real, and organize based on
#drug | real | permutation 1-1000
# towrite=tibble()
# for(disease in unique(clustermatrix$Disease)){
#   #Can filter down to only drugs I care about, aka the ones hit by the real biogrid
#   f = filter(drugscores,Disease==disease) %>% select(-Disease)
#   drugsicareabout=(f %>% filter(Network=="biogrid"))$Drugs
#   f=f %>% filter(Drugs %in% drugsicareabout) %>% arrange(Network,Drugs)
#   if(nrow(f)==0){
#     next
#   }
#   #before I pivot, f has drugs only found by real. So I should replace the latter 3
#   #cols in f with values from network
#   ftopivot = f %>% select(Network,Drugs,Scores) %>% arrange(Drugs)
#   pivoted=pivot_wider(ftopivot,values_from=Scores,names_from=Network) %>% replace(is.na(.),0) %>% arrange(Drugs)
#   
#   
#   print(head(colnames(pivoted),2))
#   #no drugs found for real, so add it in
#   #if(colnames(pivoted)[2]!="biogrid"){
#   #  pivoted$biogrid=0
#   #  pivoted=relocate(pivoted,biogrid,.after=Drugs)
#   #}
#   print(head(colnames(pivoted),2))
#   #Finding columns that have greater value than what is in col 2
#   rows=1:nrow(pivoted)
#   pvals=mclapply(rows,getpvalue,mc.cores=detectCores()-1)
#   pivoted$pvals=as.double(unlist(pvals))
#   pivoted$FDR=p.adjust(pivoted$pvals,method="BH")
#   pivoted=pivoted %>% relocate(pvals,.after=Drugs) %>% relocate(FDR,.after=pvals)
#   print(disease)
#   print(pivoted)
#   print(summary(pivoted$pvals))
#   print(summary(pivoted$FDR))
#   pivoted= pivoted %>%left_join((f %>% filter(Network=="biogrid") %>% select(Drugs,EntrezFound,SymbolFound,Clustersforgenes)),by="Drugs") %>%
#     relocate(EntrezFound,.after=Drugs) %>% relocate(SymbolFound,.after=EntrezFound) %>% relocate(Clustersforgenes,.after=SymbolFound)
#   toadd = pivoted %>% select(Drugs,EntrezFound,SymbolFound,Clustersforgenes,pvals,FDR,biogrid)
#   if(!(disease %in% c("Faketrait_2","Faketrait_5","Faketrait_7","Height_gene_shot","Left_Hand_Phone","Right_Hand_Phone"))){
#     towrite=rbind(towrite,tibble(Disease=disease,toadd))
#   }
#   #write_tsv(pivoted,paste0("/mnt/research/compbio/krishnanlab/projects/chronic_inflammation/drugs/output/drugscores/",disease,"_drugs.tsv"))
# }

#loop for getting scores if I use 1000 fake traits instead
towrite=tibble()
for(disease in unique(clustermatrix$Disease)){
  f=filter(realtraitscores,Trait==disease)
  if(nrow(f)==0){
    next
  }
  drugsicareabout=f$Drugs
  faketraitdrugs=filter(faketraitscores, Drugs %in% drugsicareabout)
  ftopivot = f %>% select(Trait,Drugs,Scores) %>% arrange(Drugs)
  faketraitdrugs = faketraitdrugs %>% select(Trait,Drugs,Scores) %>% arrange(Drugs)
  #get faketraitdrugs that match with real disease
  faketraitdrugs = faketraitdrugs %>% filter(grepl(disease,Trait))
  
  topivot=rbind(ftopivot,faketraitdrugs)
  pivoted=pivot_wider(topivot,values_from=Scores,names_from=Trait) %>% replace(is.na(.),0) %>% arrange(Drugs)
  rows=1:nrow(pivoted)
  pvals=mclapply(rows,getpvalues,5000,mc.cores=detectCores()-1)
  pivoted$pvals=as.double(unlist(pvals))
  pivoted$FDR=p.adjust(pivoted$pvals,method="BH")
  pivoted=pivoted %>% relocate(pvals,.after=Drugs) %>% relocate(FDR,.after=pvals)
  print(disease)
  print(pivoted)
  print(summary(pivoted$pvals))
  print(summary(pivoted$FDR))
  
  pivoted= pivoted %>%left_join((f %>% select(Drugs,EntrezFound,SymbolFound,Clustersforgenes)),by="Drugs") %>%
    relocate(EntrezFound,.after=Drugs) %>% relocate(SymbolFound,.after=EntrezFound) %>% relocate(Clustersforgenes,.after=SymbolFound)
  colnames(pivoted)[7]="Score"
  toadd = pivoted %>% select(Drugs,EntrezFound,SymbolFound,Clustersforgenes,pvals,FDR,Score)
  
  #if(!(disease %in% c("Faketrait_2","Faketrait_5","Faketrait_7","Height_gene_shot","Left_Hand_Phone","Right_Hand_Phone"))){
  towrite=rbind(towrite,tibble(Disease=disease,toadd))
  #}
}




#GEtting indications
dc=select(drugcentral,STRUCT_ID,DRUG_NAME)  %>% distinct
colnames(dc)[2]="Drugs"
towrite2=towrite %>% left_join(dc,by="Drugs")
ind=read_csv("../data/drugcentral/databasecsvs/indications.csv") %>% filter(relationship_name=="indication" | relationship_name=="off-label use") %>% select(-id,-concept_id)
colnames(ind)[1]="STRUCT_ID"
ind = ind %>% filter(STRUCT_ID %in% towrite2$STRUCT_ID)

indicationst=tibble()
for(i in unique(towrite2$STRUCT_ID)){
  f=filter(ind,STRUCT_ID==i)
  if(nrow(f)>0){
    indicationst=rbind(indicationst,tibble(STRUCT_ID=i,indications=paste0(f$concept_name,collapse="|"),relationship=paste0(f$relationship_name,collapse="|")))
  }else{
    indicationst=rbind(indicationst,tibble(STRUCT_ID=i,indications=NA,relationship=NA))
  }
}
disinds=towrite2 %>% left_join(indicationst,by="STRUCT_ID") #%>% select(-STRUCT_ID)

#Get gene 


#disinds=towrite2 %>% left_join(ind,by="STRUCT_ID") %>% select(-snomed_full_name,-cui_semantic_type,-snomed_conceptid,-umls_cui)
write_tsv(disinds,"/mnt/research/compbio/krishnanlab/projects/chronic_inflammation/results/GenePlexus_Drugs.tsv")

#write_tsv(networkdata,"/mnt/research/compbio/krishnanlab/projects/chronic_inflammation/drugs/output/finalnetworkdata.tsv")
#drugcentral expert curated number of times each drug appears histogram
#hist(table(drugcentral$DRUG_NAME),breaks=500)