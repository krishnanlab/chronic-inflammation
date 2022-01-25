#Order to run things in and will be a broad pipeline for how to run things

#Data setup
#prep_disease_gene_dfs.R

#prep_edgelist_hotnet2.R

#For real data, makeNodeNeighborhoodList.R and expandAndCluster.R
makeNodeNeighborhoodlist(edgelist="",directed=F,permuted=F,network="biogrid",outdir="")


#Make permutations


#for fake data, makeNodeNeighborhoodList.R and expandAndCluster.R
makeNodeNeighborhoodlist(edgelist="",directed=F,permuted=T,network="biogrid",outdir="")
