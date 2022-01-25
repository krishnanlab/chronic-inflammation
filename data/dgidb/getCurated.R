#Call API to get expert curated interactions for each gene

library(tidyverse)
#json library
library(httr)
library(jsonlite)

#load interactions
interactions=read_tsv("interactions.tsv",col_names=T)

#get all unique gene ids from this
geneids=unique(interactions$gene_name)

#function for api call
#curl https://dgidb.org/api/v2/interactions.json?genes=TRPC6
#q=GET("https://dgidb.org/api/v2/interactions.json?genes=TRPC6&source_trust_levels=Expert%20curated")
#rawToChar(q$content)
#idnum is the iterator for geneids
getjson_tibble<-function(idnum){
  print(idnum)
  url=paste0("https://dgidb.org/api/v2/interactions.json?genes=",geneids[idnum],"&source_trust_levels=Expert%20curated")
  stuff=GET(url)
  #in a json string, I have to remove the first 17 characters, and the last 42
  #This pattern is always the same as for this use case I will always find terms, meaning
  #never ambiguous or unmatched
  r=rawToChar(stuff$content)
  endgone=substring(r,1,nchar(r)-42)
  begingone=substring(endgone,18,nchar(endgone))
  lastcharacter=substring(begingone,nchar(begingone),nchar(begingone))
  #Beginning will always be the same. If ever anything wrong here, then do soemthing
  #where I find the pattern I want, then remove everything after that!
  if(lastcharacter!="}"){
    print("pattern broken, something in ambiguous or unmateched")
    cat("bad for idnum ",idnum,"Which is gene ",geneids[idnum],"\n")
  }
  finaljson=begingone
  #fromJSON(finaljson)$interactions will have what I need
  jsontib=fromJSON(finaljson)$interactions %>% as_tibble
  #Only need drugName, score, and indication on what gene this is. I can get
  #sources, pmids, and interactionTypes from interactions.tsv
  #check to see if actually have interactions
  if(nrow(jsontib)>0){
    fintib=select(jsontib,drugName,score)
    fintib$gene_name=geneids[idnum]
  }
  else{
    fintib=tibble()
  }
  return(fintib)
}

#idnum will be a number referring to i in geneids[i]
count=1
for(idnum in 1:length(geneids)){
  tib=getjson_tibble(idnum)
  #add this tib to final tib
  if(count==1){
    finaltib=tib
  }
  else{
    finaltib=rbind(finaltib,tib)
  }
  count=count+1
}
finaltib$expert="yes"
write_tsv(finaltib,"expertcuratedinteractions.tsv",col_names=T)

#fromJSON(jsons[1],simplifyVector=F)
#w=fromJSON(jsons[1])$matchedTerms
#toJSON(w) %>%spread_all

#got rid of the matched term crap and it loads much better
#fromJSON("temp.json")



