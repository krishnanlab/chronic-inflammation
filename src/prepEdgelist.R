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
#' @args[4] weights are in V3 and I want to keep them T/F

library(tidyverse)
library(igraph)
source("../src/chronic_inflammation_functions.R")
args <- commandArgs(TRUE)

# load edgelist
edge = read.delim(args[1], header = F)
out_dir  = args[2]
net_name = args[3]

# make index file
V1 = unique(edge$V1)
V2 = unique(edge$V2)

ids = unique(c(V1, V2))
ids_human = entrez2symbol(as.character(ids))

ids = as.data.frame(ids_human)
colnames(ids) = "SYMBOL"
ids$ENTREZ = names(ids_human)
ids$index = 1:nrow(ids)

ids = 
  ids %>%
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
            file = paste0(out_dir, "/", net_name, "_human_entrez_and_gene_symbols.txt"),
            sep = "\t",
            row.names = F,
            col.names = F,
            quote = F)

# filter non-human genes from
# edge list
edge_human = 
  edge %>%
  filter(V1 %in% ids$ENTREZ &
           V2 %in% ids$ENTREZ)

write.table(edge_human, 
            file = paste0(out_dir, "/", net_name, "_human_entrez_edgelist.txt"),
            sep = "\t",
            row.names = F,
            col.names = F,
            quote = F)

# make igraph object
edge_human$V1 = as.character(edge_human$V1)
edge_human$V2 = as.character(edge_human$V2)

edge_mat = as.matrix(edge_human[c(1,2)])

g = graph_from_edgelist(edge_mat, directed = FALSE)

if(as.logical(args[4]) == TRUE){
  # add weight
  E(g)$weight = edge_human$V3
}

#save
save(g, file = paste0(out_dir, "/", net_name, "_igraph.Rdata"))
  
print("all done!")

