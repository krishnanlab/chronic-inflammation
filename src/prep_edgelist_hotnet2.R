# prep index and edgelist files 
# the edgelist needed for hotnet2 uses indices, not gene names
# https://github.com/raphael-group/hotnet2/blob/master/example/example_edgelist.txt
# the edgelist I have for biogrid is entrez IDs
# make an index to gene file 
# filter out non-human genes
# https://github.com/raphael-group/hotnet2/blob/master/example/example_gene_index.txt
# then replace the gene names in the edgelist with the indices
#' @args[1] path to tab delimited edgelist
#' @args[2] path to output dir
#' @args[3] network name
#' @args[4] make index edgelist T/F

library(tidyverse)
source("chronic_inflammation_functions.R")
args <- commandArgs(TRUE)

# load edgelist
edge = read.delim(args[1], header = F)
out_dir  = args[2]
net_name = args[3]
args[4]=F
# make index file
V1 = unique(edge$V1)
V2 = unique(edge$V2)

ids = unique(c(V1, V2))
ids_human = entrez2symbol(as.character(ids))

ids = as.data.frame(ids_human)
colnames(ids) = "SYMBOL"
ids$ENTREZ = names(ids_human)
ids$index = 1:nrow(ids)

ids = ids %>%
  dplyr::select(index, 
                ENTREZ, 
                SYMBOL)
  
write.table(ids %>%
              dplyr::select(index,
                            ENTREZ), 
            file = paste0(out_dir, "/", net_name, "_gene_index.txt"),
            sep = "\t",
            row.names = F,
            col.names = F,
            quote = F)

write.table(ids,
            file = paste0(out_dir, "/", net_name, "_index_with_gene_symbols.txt"),
            sep = "\t",
            row.names = F,
            col.names = F,
            quote = F)

# filter non-human genes from
# edge list
edge_human = edge %>%
  filter(V1 %in% ids$ENTREZ &
           V2 %in% ids$ENTREZ)

write.table(edge_human, 
            file = paste0(out_dir, "/", net_name, "_entrez_edgelist.txt"),
            sep = "\t",
            row.names = F,
            col.names = F,
            quote = F)

if(args[4] == TRUE){# replace entrex ids with index incase hotnet2 needs this
  edge_index = edge_human
  edge_index$V1 = paste0("id", edge_index$V1)
  edge_index$V2 = paste0("id", edge_index$V2)
  
  ids$ENTREZ = paste0("id", ids$ENTREZ)
  
  for(i in 1:nrow(ids)){
    
    print(i)
    
    edge_index$V1 = gsub(paste0('\\<', ids[i,2], '\\>'),
                         ids[i,2],
                         edge_index$V1)
    
    edge_index$V2 = gsub(paste0('\\<', ids[i,2], '\\>'),
                         ids[i,1],
                         edge_index$V2)
  }
  
  print("index edgelist done")
  
  write.table(edge_index, 
              file = paste0(out_dir, "/", net_name, "_index_edgelist.txt"),
              sep = "\t",
              row.names = F,
              col.names = F,
              quote = F)
}

print("all done!")

