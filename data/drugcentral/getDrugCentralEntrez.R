#Here I get drugcenral entrez which I will use for drug searching. I also only use FDA drugs as well
#as getting ATC codes

library(tidyverse)

#Now reading the files I got frm the database
#acttable has targets
targets=read_csv("./databasecsvs/acttable.csv")

#Getting only fda approved from approval
approval=read_csv("./databasecsvs/approval.csv")

#Use only preferred names
ids=read_csv("./databasecsvs/synonyms.csv") %>% filter(id %in% approval$struct_id) %>% filter(preferred_name==1)
#atc code. FDA drugs do not necessarily have an ATC code, while many have multiple
atcs=read_csv("./databasecsvs/struct2atc.csv") %>% dplyr::select(-id)

#subset targets for if gene not na, and homo sapien
targets= targets %>% filter(organism=="Homo sapiens") %>% filter(!is.na(gene)) %>% dplyr::select(struct_id,target_class,gene,tdl)

colnames(ids)[2]="struct_id"
ids=ids %>% dplyr::select(struct_id,name)
#collapse name column in ids for each struct id
colids=tibble()
for(i in unique(ids$struct_id)){
  f=filter(ids,struct_id==i)
  colids=rbind(colids,tibble(struct_id=i,names=paste0(f$name,collapse=", ")))
}
#get only FDA approved
colids=colids %>% filter(struct_id %in% approval$struct_id) %>% arrange(struct_id)
#collape atcs
colatcs=tibble()
for(i in unique(ids$struct_id)){
  f=filter(atcs,struct_id==i)
  colatcs=rbind(colatcs,tibble(struct_id=i,atcs=paste0(f$atc_code,collapse=", ")))
}
#join with ATC code
colids=colids %>% left_join(colatcs,by="struct_id")
#join with targets
colids=colids %>% left_join(targets,by="struct_id")
colnames(colids)=c("STRUCT_ID","DRUG_NAME","ATC_CODE","TARGET_CLASS","GENE","TDL")
drugcentral=colids %>% filter(!is.na(GENE))
#drugcentral=read_tsv("../data/drugcentral/drug.target.interaction.tsv")
#drugcentral=drugcentral %>% filter(ORGANISM=="Homo sapiens")

library(org.Hs.eg.db)
entrezvector=unlist(as.list(org.Hs.egSYMBOL2EG))
symentrez=as_tibble(entrezvector)
symentrez$GENE=names(entrezvector)
colnames(symentrez)[1]="Entrez"
#At this point, need to deal with messiness. Capitalize all gene symbol characters
drugcentral$GENE=toupper(drugcentral$GENE)
drugcentral=drugcentral %>% separate_rows(GENE,TDL,sep="\\|")

#now join
drugcentral2=left_join(drugcentral,symentrez,by="GENE") %>% filter(!is.na(Entrez)) %>% filter(TDL %in% c("Tclin","Tchem"))%>%
  dplyr::select(STRUCT_ID,DRUG_NAME,GENE,ATC_CODE,TARGET_CLASS,TDL,Entrez)
drugcentral2[drugcentral2==""]<-NA

write_tsv(drugcentral2,"drugcentral_entrez.tsv",col_names=T)

#17 genes don't have an entrez
table(filter(drugcentral2,is.na(Entrez))$GENE)
