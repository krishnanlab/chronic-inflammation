library(tidyverse)
expertcurated=read_tsv("expertcuratedinteractions.tsv")
#expertcurated=expertcurated %>% filter(ORGANISM=="Homo sapiens")
#expertcurated= dplyr::select(expertcurated,DRUG_NAME,GENE)

library(org.Hs.eg.db)
entrezvector=unlist(as.list(org.Hs.egSYMBOL2EG))
symentrez=as_tibble(entrezvector)
symentrez$gene_name=names(entrezvector)
colnames(symentrez)[1]="Entrez"
#At this point, need to deal with messiness. Capitalize all gene symbol characters
expertcurated$gene_name=toupper(expertcurated$gene_name)
#Second thing is that you will have Gene1 | Gene2 | Gene3 ... | Gene N. How deal?
#Split on | and create new lines like I do with loading other conversions
#testframe=data.frame(drugName=rep(expertcurated$drugName,sapply(split,length)),gene_name=unlist(split))
#testframe %>% mutate_if(is.factor,as.character) %>% as_tibble(.)-> expertcurated2

#now join
expertcurated2=left_join(expertcurated,symentrez,by="gene_name")
write_tsv(expertcurated2,"expertcurated_entrez.tsv",col_names=T)

#17 genes don't have an entrez
table(filter(expertcurated2,is.na(Entrez))$Entrez)
