# loadRData ---------------------------------------------------------------
#' @name loadRData
#' @description loads an RData file, and returns it
#' @param fileName path to RData file
# https://stackoverflow.com/questions/5577221/how-can-i-load-an-object-into-a-variable-name-that-i-specify-from-an-r-data-file

loadRData <- function(fileName){
  #loads an RData file, and returns it
  load(fileName)
  get(ls()[ls() != "fileName"])
}

# makeUndirectedEdgelist --------------------------------------------------
#' @name makeUndirectedEdgelist
#' @description Takes an edgelist df, switches the to and from node
#' columns and adds those rows back to the original edgelist
#' @param edgelist df with two col
#' @return new edgelist where every "from" node is now also 
#' a "to" node and vice versa

makeUndirectedEdgelist = function(edgelist){
  switch = edgelist[c(2,1)]
  colnames(switch) = colnames(edgelist)
  new_edgelist = rbind(edgelist, switch)
  new_edgelist = unique(new_edgelist)
  return(new_edgelist)
}

# makeNodeNeighborhoodList ------------------------------------------------
#' @name makeNodeNeighborhoodList
#' @description Takes an edgelist df and returns a list of 
#' first degree neighbors for every node
#' @param edgelist df with two col
#' @param directed TRUE/FALSE 
#' @return named list of first degree neighbors for every node

makeNodeNeighborhoodList = function(edgelist, directed = FALSE){
  
  require(dplyr) 
  
  # if undirected flip the cols and rbind
  if(directed == FALSE){
    nh = makeUndirectedEdgelist(edgelist)
  } else {
    nh = edgelist
  }
  
  colnames(nh) = c("Gene", "DEFINITION")
  nh_def = unique(nh$DEFINITION)
  ln = length(nh_def)
  nh_list = list()
  
  print("finding node neighborhoods")
  for (i in 1:ln){
    nh_list[[i]] = nh %>%
      filter(DEFINITION == nh_def[i]) %>%
      pull(Gene)
  }
  names(nh_list) = nh_def
  return(nh_list)
}


# expandGeneList ----------------------------------------------------------
#' @name expandGeneList
#' @description performs a hypergeometric test between genes of interest 
#' and the first degree neighbors of every gene in a network.
#' If a node's neighbors are enriched for the genes of interest, the node 
#' becomes a new gene of interest. Process is repeated for a specifed number
#' of iterations.
#' @param goi_df a data frame with the genes of interest in column 1 and 
#' the associated definition (i.e. disease) in column 2
#' @param neighbor_list the output of makeNodeNeighborhoodList
#' a list of first degree neighbors of each node in a network
#' @param cutoff the pval cutoff for enrichment
#' @param FDR_TF TRUE or FALSE weither to use FDR. 
#' FALSE used uncorrected pval
#' @param iterations number of times to repeat the overlap test
#' @return returns goi_df but with new genes of interest added

expandGeneList <- function(goi_df, 
                           neighbor_list, 
                           cutoff, 
                           FDR_TF,
                           iterations){
  
  require(dplyr)
  
  net_genes = names(neighbor_list)
  bg = length(net_genes)
  
  goi_df = goi_df %>%
    filter(Gene %in% net_genes)
  
  nGroup2 = unlist(lapply(neighbor_list, length))
  
  for(i in 1:iterations){
  
    print(paste0("iteration ", i))
    colnames(goi_df) = c("Gene","DEFINITION")
    nGroup1 = nrow(goi_df)
    INT = lapply(neighbor_list, function(x) x[x %in% goi_df$Gene])
    nINT = unlist(lapply(INT, length))
  
    df = cbind(nGroup2, nINT)
  
    print("testing neighborhoods for enrichment")
    
    apply_phyper = function(df_row, background_n, nGroup1){
      phyper(df_row[2]-1, 
            df_row[1], 
            background_n-df_row[2], 
            nGroup1,
            lower.tail= FALSE)}
  
    pvals = apply(df, 1, apply_phyper, bg, nGroup1)
    FDR = p.adjust(pvals, method = "BH")
  
    if(FDR_TF==TRUE){
      new_nodes = names(FDR)[which(FDR < cutoff)]
    } else {
      new_nodes = names(pvals)[which(pvals < cutoff)]}
  
    add_nodes = new_nodes[!new_nodes %in% goi_df$Gene]
    print(paste0("n new nodes = ", length(add_nodes)))
    if(length(add_nodes) < 1) {break}
    add_nodes = as.data.frame(add_nodes)
    add_nodes$DEFINITION = unique(goi_df$DEFINITION)
    colnames(add_nodes) = colnames(goi_df)
    goi_df = rbind(goi_df, add_nodes)
    last_i = i
  }
  out = list(goi_df, last_i)
  names(out) = c("genes", "n_iterations")
  return(out)
}

# findEnrichedGOBP --------------------------------------------------------
#' @name findEnrichedGOBP
#' @description uses topGO to find enriched GOBP terms for a list of genes.
#' ------------ Pvalues from weight01 fisher test
#' @param goi character vector of genes of interest
#' @param background character vector of all background genes
#' @param org.db the org.db for the organism of interest from bioconductor. 
#' ------------- Must already be installed.
#' @param id_type gene identifier to use. 
#' ------------- c("entrez","genbank","alias","ensembl","symbol","genename")
#' @param min_size smallest number of genes for a GO term to be used
#' @param max_size largest number of genes for a GO term to be used
#' @return a data frame with the pvalue and overlapping genes 
#' ------- between the genes of interest and GOBPs. 

findEnrichedGOBP = function(goi, 
                            background, 
                            org.db, 
                            id_type,
                            min_size,
                            max_size){
  
  require(paste(org.db), character.only = TRUE)
  require(tidyverse)
  require(topGO)
  
  submit = factor(ifelse(background %in% goi, 1, 0))
  names(submit) = background
  
  GO = new("topGOdata",
           ontology = "BP",
           allGenes = submit,
           nodeSize = min_size,
           annot = annFUN.org,
           mapping = org.db,
           ID = id_type)
  
  nTerms = length(usedGO(GO))
  
  Fisher <- runTest(GO,
                    algorithm = "weight01", 
                    statistic = "fisher")
  
  results = GenTable(GO, 
                     pval = Fisher, 
                     topNodes = nTerms,
                     numChar=1000)
  
  results = results %>%
    filter(Annotated <=  max_size)
  
  results$FDR = p.adjust(results$pval, method = "BH")
  
  ann.genes = genesInTerm(GO, results$GO.ID)
  symbol = names(submit)[submit == 1]
  
  results$OverlappingGenes = sapply(results$GO.ID, 
                                    function(x){
                                    sig = ann.genes[[x]][ann.genes[[x]] %in% symbol]
                                    sig = paste(sig, collapse = ", ")
                                    return(sig)}
                                    )
  return(results)
}


# overlapSets ------------------------------------------------------------
#' @name overlap_sets
#' @description performs pairwise hypergeometric tests between mutiple sets 
#' of genes i.e. cluster marker genes from two different single cell
#' different single cell data sets
#' @param table1 a dataframe with two columns, 
#' ie gene and cluster assignment, from dataset 1      
#' @param table2 a dataframe with two columns, 
#' ie gene and cluster assignment, from dataset 2
#' format example
#' --------------  
#' Gene  Cluster
#' Gene1 cluster1
#' Gene2 cluster1
#' Gene3 cluster2
#' Gene4 cluster3
#' @param background the number of genes in the gene universe 
#' ie. all expressed genes
#' @return a list with two elements: 
#' [1]intersecting.genes: a tidy table where each row shows two groups 
#' (one from table1 and one from table2) and an intersecting gene or "none"
#' [2] pval.scores: a tidy table where each row is a pair of groups 
#' (one from table1 and one from table2) with the number of genes in each group,
#' the nunber of intersecting genes, the uncorrected pval, and the enrichment score

overlapSets <- function(table1, table2, background){
  require(dplyr)
  
  # make table2 into a list divided by group
  colnames(table1) = c("Gene","Group")
  tab1_groups = unique(table1$Group)
  tab1_list = list()
  
  for(group in tab1_groups){
    tab1_list[[group]] = table1 %>%
      filter(Group == group)
  }
  
  # make table2 into a list divided by group
  colnames(table2) = c("Gene","Group")
  tab2_groups = unique(table2$Group)
  tab2_list = list()
  
  for(group in tab2_groups){
    tab2_list[[group]] = table2 %>%
      filter(Group == group)
  }
  
  # make a list of dfs with the names of genes shared between
  # every group in table1 and every group in table2
  INT = list()
  # make a list of dfs with the number of genes shared between
  # every group in table1 and every group in table2
  nINT = list()
  
  for(group in tab1_groups){
   genes = tab1_list[[group]]$Gene
   INT_list = lapply(tab2_list, 
                         function(x){
                           y = intersect(x$Gene, genes);
                           if(length(y) > 0){
                           z = data.frame(SharedGene = y)}
                           else{z = data.frame(SharedGene = "none")};
                           z
                         })
   
   for(i in 1:length(INT_list)){
     INT_list[[i]]$Group2 = names(INT_list)[i]
   }

   INT[[group]] = do.call(rbind,INT_list)
   INT[[group]]$Group1 = group
   
   INT[[group]] = INT[[group]] %>%
     select(Group1, Group2, SharedGene)
   
   nINT_list = lapply(tab2_list, 
                     function(x){
                       y = intersect(x$Gene, genes);
                       z = length(y);
                       z
                     })
   
   nINT_vec = unlist(nINT_list)
   nGroup2 = unlist(lapply(tab2_list, nrow))
   
   nINT[[group]] = data.frame(Group1 = group,
                              nGroup1 = length(genes),
                              Group2 = names(nINT_vec),
                              nGroup2  = nGroup2,
                              nSharedGenes = nINT_vec)
   
   
  }
  
  INTdf = do.call(rbind, INT)
  nINTdf = do.call(rbind, nINT)
  
  nINTdf = nINTdf %>%
    mutate(Pval = phyper(nSharedGenes-1,
                         nGroup2,
                         background-nGroup2,
                         nGroup1,
                         lower.tail= FALSE)
                         )
  
  nINTdf = nINTdf %>%
    mutate(Enrichment = 
             log2(
               (nSharedGenes/nGroup2)/
                    (nGroup1/background))
           )
  
  return(list(intersecting.genes = INTdf, 
              pval.scores = nINTdf))
}  

# entrez2symbol -----------------------------------------------------------
entrez2symbol = function(entrez_ids){
  require("org.Hs.eg.db")
  x <- org.Hs.egSYMBOL
  included_entrez <- names(as.list(x))
  included_oi = entrez_ids[entrez_ids %in% included_entrez]
  xx <- as.list(x[included_oi])
  gene_symbols = unlist(xx)
  return(gene_symbols)
  unloadNamespace("org.Hs.eg.db")
  detach("package:AnnotationDbi")
}

# symbol2entrez -----------------------------------------------------------
symbol2entrez = function(gene_symbols){
  require("org.Hs.eg.db")
  x <- org.Hs.egSYMBOL2EG
  included_symbols <- names(as.list(x))
  included_oi = gene_symbols[gene_symbols %in% included_symbols]
  xx <- as.list(x[included_oi])
  entrez_ids = unlist(xx)
  return(entrez_ids)
  unloadNamespace("org.Hs.eg.db")
  detach("package:AnnotationDbi")
}