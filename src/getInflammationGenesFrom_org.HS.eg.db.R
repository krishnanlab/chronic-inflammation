# programatically access genes associated with gene ontology terms
library(org.Hs.eg.db) #org.Hs.eg.db_3.13.0
disease_gene_dir = "../data/disease_gene_files"

# experimental evidence codes
# http://geneontology.org/docs/guide-go-evidence-codes/
expr = c("EXP", #Inferred from Experiment (EXP)
         "IDA", #Inferred from Direct Assay (IDA)
         "IPI", #Inferred from Physical Interaction (IPI)
         "IMP", #Inferred from Mutant Phenotype (IMP)
         "IGI", #Inferred from Genetic Interaction (IGI)
         "IEP") #Inferred from Expression Pattern (IEP)

# high throughput experiments
ht_expr = c("HTP", #Inferred from High Throughput Experiment (HTP)
            "HDA", #Inferred from High Throughput Direct Assay (HDA)
            "HMP", #Inferred from High Throughput Mutant Phenotype (HMP)
            "HGI", #Inferred from High Throughput Genetic Interaction (HGI)
            "HEP") #Inferred from High Throughput Expression Pattern (HEP)

# chronic inflammation no propagation from child terms
GO2EG <- as.list(org.Hs.egGO2EG)
go_genes <- GO2EG["GO:0002544"][[1]] #chronic inflammatory response

df = as.data.frame(unique(go_genes))
colnames(df) = "Gene"
df$Disease = "chronic_inflammatory_response_GO2EG"
write.table(df, 
            file = paste0(disease_gene_dir, "/chronic_inflammatory_response_GO2EG.txt"),
            row.names = F,
            quote = F,
            sep = "\t")

# chronic inflammation with propagation from child terms
GO2ALLEGS = as.list(org.Hs.egGO2ALLEGS)
go_all_genes = GO2ALLEGS["GO:0002544"][[1]] #chronic inflammatory response

df = as.data.frame(unique(go_all_genes))
colnames(df) = "Gene"
df$Disease = "chronic_inflammatory_response_GO2ALLEGS"
write.table(df, 
            file = paste0(disease_gene_dir, "/chronic_inflammatory_response_GO2ALLEGS.txt"),
            row.names = F,
            quote = F,
            sep = "\t")

# inflammation no propagation from child terms
go_genes <- GO2EG["GO:0006954"][[1]] #inflammatory response

#filter for experimental evidence
go_genes = go_genes[which(names(go_genes) %in% c(expr, ht_expr))]

df = as.data.frame(unique(go_genes))
colnames(df) = "Gene"
df$Disease = "inflammatory_response_GO2EG_expr"
write.table(df, 
            file = paste0(disease_gene_dir, "/inflammatory_response_GO2EG_expr.txt"),
            row.names = F,
            quote = F,
            sep = "\t")

# inflammation with propagation from child terms
go_all_genes = GO2ALLEGS["GO:0006954"][[1]] #inflammatory response

#filter for experimental evidence
go_all_genes = go_all_genes[which(names(go_all_genes) %in% c(expr, ht_expr))]

df = as.data.frame(unique(go_all_genes))
colnames(df) = "Gene"
df$Disease = "inflammatory_response_GO2ALLEGS_expr"
write.table(df, 
            file = paste0(disease_gene_dir, "/inflammatory_response_GO2ALLEGS_expr.txt"),
            row.names = F,
            quote = F,
            sep = "\t")
