library(tidyverse)
setwd("/mnt/ufs18/rs-027/compbio/krishnanlab/data/disease-gene_annotations/disgenet/2022_04")

# map between UMI ids and other ids for synonmyns 
disease_map = read.delim("disease_mappings.tsv.gz")

# find diseases with at least one associated gene
gene2disease = read.delim("all_gene_disease_associations.tsv.gz")

# filter out diseaseType == "phenotype"
# filter for diseaseSemanticType of interest
# filter out assoications based only on text mining the literature
# https://www.disgenet.org/dbinfo
lit_only = c("LHGDN", "BEFREE", "BEFREE;LHGDN")
# filter for diseaseSemanticType of interest
keep_diseaseSemanticType = c("Mental or Behavioral Dysfunction",
                             "Disease or Syndrome",
                             "Neoplastic Process")

gene2disease_filt = 
  gene2disease %>%
  filter(!diseaseType == "phenotype", # filter out diseaseType == "phenotype"
         !source %in% lit_only,
         diseaseSemanticType %in% keep_diseaseSemanticType)

# pull disease ids from gene2disease_filt
gene2disease_ids = 
  gene2disease_filt %>%
  select(diseaseId) %>%
  distinct() %>%
  pull(diseaseId)

length(gene2disease_ids) #11093
sum(gene2disease_ids %in% disease_map$diseaseId) #10911

# get ids missing from the mapping file
gene2disease_ids_not_mapped = gene2disease_ids[!gene2disease_ids %in% disease_map$diseaseId]

# get disease names of unmapped ids
gene2disease_ids_not_mapped_names = 
  gene2disease_filt %>%
  filter(diseaseId %in% gene2disease_ids_not_mapped) %>%
  select(diseaseId, diseaseName) %>%
  distinct()

# add vocabulary, code, vocabularyName cols
# filled with NA so we can bind these with the 
# disease_map table later
colnames(gene2disease_ids_not_mapped_names) = c("diseaseId", "name")

gene2disease_ids_not_mapped_names = 
  gene2disease_ids_not_mapped_names %>%
  mutate(vocabulary = NA,
         code = NA,
         vocabularyName = NA)
  
# find diseases with at least one associated drug
indications = read.delim("/mnt/research/compbio/krishnanlab/projects/chronic_inflammation/drugs/data/drugcentral/indications.tsv")

# filter for indication or off-label use
# filter for diseases and syndromes T047, Neoplastic process T191, Mental or behavioral dysfunction T048
# https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4450611/
indications_filt = 
  indications %>%
  filter(relationship_name %in% c("indication","off-label use"),
         cui_semantic_type %in% c("T047", "T191", "T048"))


# pull disease ids from indications_filt
indications_id = 
  indications_filt %>%
  select(umls_cui) %>%
  distinct() %>%
  pull(umls_cui)
  
length(indications_id) # 1310

# how many indications have synonyms in disease_map
sum(indications_id %in% disease_map$diseaseId) # 961

# get ids missing from the mapping file
indication_ids_not_mapped = indications_id[!indications_id %in% disease_map$diseaseId]

# use "snomed" as vocabulary
# use the snomed_conceptid concept id as code
# use snomed_full_name as vocabularyName

indication_ids_not_mapped_names = 
  indications_filt %>%
  filter(umls_cui %in% indication_ids_not_mapped) %>%
  mutate(vocabulary = "snomed") %>%
  select(umls_cui, 
         concept_name, 
         vocabulary,
         snomed_conceptid, 
         snomed_full_name) %>%
  distinct() 

colnames(indication_ids_not_mapped_names) = colnames(disease_map)

# remove gene2disease_ids_not_mapped that 
# are also in indication_ids_not_mapped
# from gene2disease_ids_not_mapped_names

gene2disease_ids_not_mapped_names = 
  gene2disease_ids_not_mapped_names %>%
  filter(!diseaseId %in% indication_ids_not_mapped)

# filter disease_map for cuis of intersest
id_oi = unique(c(gene2disease_ids, indications_id))

disease_map_filt = 
  disease_map %>%
  filter(diseaseId %in% id_oi)
  
# add ids not in mapping file
disease_map_filt = rbind(disease_map_filt, indication_ids_not_mapped_names)
disease_map_filt = rbind(disease_map_filt, gene2disease_ids_not_mapped_names)

# save for Matt
write.csv(disease_map_filt, file = "/mnt/research/compbio/krishnanlab/projects/chronic_inflammation/data/diseases_with_synonmyns_for_matt.csv")



