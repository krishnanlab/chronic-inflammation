#' @param args[1] path to edgelist (no header)
#' @param args[2] directed T/F
#' @param args[3] permuted T/F
#' @param args[4] netname
#' @param args[5] outdir

library(dplyr)
source("chronic_inflammation_functions.R")
args <- commandArgs(TRUE)

outdir = args[5]
if (!dir.exists(outdir)) {dir.create(outdir)}
#setwd(outdir)

if(args[3] == TRUE){
  edgelist = read.delim(args[1], header = F, sep = " ")
  edgelist$V3 = NULL
} else {edgelist = read.delim(args[1], header = F)}

nh_list = makeNodeNeighborhoodList(edgelist, as.logical(args[2]))

netname = args[4]
save(nh_list, file = paste0(outdir, 
                            "/", 
                            netname, 
                            "_neighbor_list.Rdata"))
