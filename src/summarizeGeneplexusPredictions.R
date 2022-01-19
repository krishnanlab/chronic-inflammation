#' @param args[1] path to dir with gene plexus predictions
#' @param args[2] output dir
#' @param args[3] average cv threshold

# setup -------------------------------------------------------------------
library(tidyverse)
library(grid)
library(mccf1)
args <- commandArgs(TRUE)

# load disease gene predictions
pred_files = list.files(args[1], 
                        pattern = "predictions.tsv", 
                        full.names = TRUE)

#out_dir = "/mnt/research/compbio/krishnanlab/projects/chronic_inflammation/results/GenePlexus_parameter_checks"
out_dir = args[2]

getGeneplexusSummary <- function(pred_file){
  
  # load pred file
  pred = read.delim(pred_file)
  # get file info
  file_sep = str_split(basename(pred_file), pattern = "--")
  file_sep = unlist(file_sep)
  disease = file_sep[1]
  net = file_sep[2]
  features = file_sep[3]
  GSC = file_sep[4]
  print(paste0(disease, "--", net, "--", features))
  
  # get cv file
  cv_file = gsub("predictions.tsv", "CVvalues.txt", pred_file)
  # load cv results
  if(file.exists(cv_file)){
    cv = read.delim(cv_file)
    # get average log2(auPRC/prior)
    avgP = mean(cv[1,])
  } else {avgP = NA}

  #pick prediction threshold with mccf1
  response = 
    pred %>%
    filter(!Class.Label == "U") %>%
    mutate(Response = 
             case_when(Class.Label == "P" ~ 1,
                       Class.Label == "N" ~ 0)) %>%
    pull(Response)
  
  predictor = 
    pred %>%
    filter(!Class.Label == "U") %>%
    pull(Probability)
  
  mccf1_object = mccf1(response, predictor)
  mccf1_summary = summary(mccf1_object)
  best_threshold = mccf1_summary[1,2]
  
  disease_genes = 
    pred %>%
    filter(Probability >= best_threshold) %>%
    pull(Entrez)
    
  # tally number of seed genes
  n_seeds = nrow(pred %>% filter(Class.Label == "P"))
  # tally number negative genes
  n_neg = nrow(pred %>% filter(Class.Label == "N"))
  # how many new genes?
  n_new_genes = length(disease_genes) - n_seeds

  summary_df = data.frame(Disease = disease,
                          Network = net,
                          Features = features,
                          negativesFrom = GSC,
                          nSeeds = n_seeds,
                          nNegs = n_neg,
                          nNewGenes = n_new_genes,
                          Mccf1Threshold = best_threshold,
                          CVscore = avgP)
  
  return(summary_df)
}

summary_list = lapply(pred_files, getGeneplexusSummary)
summary_df = do.call(rbind, summary_list)
save(summary_df,  file =paste0(out_dir, "/Geneplexus_summary.Rdata"))

# tally number of traits with CVscore > thresh

thresh = as.numeric(args[3])

nModels = 
  summary_df %>%
  group_by(Network, 
           Features) %>%
  tally() %>%
  ungroup() %>%
  select(n) %>%
  distinct() %>%
  pull(n)

cv2_tally = 
  summary_df %>%
  group_by(Network, 
           Features) %>%
  tally(CVscore >= thresh) %>%
  mutate(Proportion = n/nModels)

#plot bar
cv2_tally %>%
  ggplot(aes(x = Features,
             y = Proportion)) +
  geom_col(position="dodge") +
  geom_text(aes(label = n), vjust = -0.2) +
  ylab(paste0("Proportion log2(auPRC/prior) >=", thresh)) + 
  theme_classic() +
  facet_wrap(~Network, ncol = 3)

ggsave(file = paste0(out_dir, "/proportion_CVscore_greater", thresh, ".pdf"),
       width = 9,
       height =4)

# box + dot CV score
ggplot(summary_df, 
       aes(x = Features, 
           y = CVscore)) + 
  geom_jitter(shape=16, 
              colour = "#90A4AE",
              position=position_jitter(0.2)) +
  geom_boxplot(color = c("#C62828"),
               fill = "NA",
               outlier.size=0, 
               outlier.shape=NA) +
  theme_classic() +
  geom_hline(yintercept=thresh, 
             linetype="dashed", 
             color = "blue") +
  ylab("average log2(auPRC/prior)") +
  facet_wrap(~Network, 
             ncol = 3)
  
  ggsave(file = paste0(out_dir, "/average_CVscore.pdf"),
         width = 9,
         height =4)
  
# box plot of prediction scores for diseases ---------------------------
# grab disease files, not ukbb files
# print cv score above
dg = read.delim("../data/disgenet_2020-10_curated_gene_disease_associations.tsv")
disease_cuid = read.delim("../data/chronic_inflammation_diseases_non-ovlp_cuid.txt", header = F)
disease_cuid = disease_cuid$V1

disease_names = 
  dg %>%
  filter(diseaseId %in% disease_cuid) %>%
  select(diseaseName) %>%
  distinct() %>%
  pull(diseaseName)

disease_names = gsub(" ", "_", disease_names)
disease_names = gsub("([_])|[[:punct:]]", "\\1", disease_names)

disease_pred_files = grep(paste(disease_names, collapse="|"), 
                          pred_files,
                          value = T)

file_sep = str_split(basename(disease_pred_files), 
                     pattern = "--")

disease_pred_list = lapply(disease_pred_files,
                           read.delim)

for(i in 1:length(disease_pred_list)){
  
  disease_pred_list[[i]]$Disease = file_sep[[i]][1]
  disease_pred_list[[i]]$Network = file_sep[[i]][2]
  disease_pred_list[[i]]$Features = file_sep[[i]][3]
}

disease_pred_df = do.call(rbind, disease_pred_list)

# add cv score info
cv2 = 
  summary_df %>%
  select(Network,
         Features,
         CVscore)

cv2$CVscore = round(cv2$CVscore, digits = 2)

all_networks = unique(disease_pred_df$Network)
all_feaures = unique(disease_pred_df$Features)

for(n in all_networks){
  for(f in all_feaures){
    
    cv2_plot = 
      cv2 %>%
      filter(Network == n,
             Features == f)
      
    disease_pred_df %>%
      filter(Network == n,
             Features == f) %>%
    ggplot(aes(x = Class.Label, 
               y = Probability)) +
      geom_boxplot() +
      geom_text(data    = cv2_plot,
                size    = 5,
                mapping = aes(x = -Inf, y = Inf, label = CVscore),
                hjust   = 1,
                vjust   = -1) +
      facet_wrap(~Disease, 
                 nrow = 9) +
      theme_classic()
     
    
    ggsave(file = paste0(out_dir, "/", n, "_", f, "_", "prediction_boxplots.pdf"),
           width = 24,
           height = 24)
    
  }
}

# just alzheimers as an example

disease_pred_df %>%
  filter(Disease == "Alzheimers_Disease") %>%
  ggplot(aes(x = Features, 
             y = Probability,
             fill = Class.Label)) +
  geom_boxplot() +
  facet_wrap(~Network, 
             nrow = 3) +
  theme_classic()

ggsave(file = paste0(out_dir, "/alzheimers_prediction_boxplots.pdf"),
       width = 9,
       height = 9)

# plot proportion new genes added for string

string = 
  summary_df %>% 
  filter(CVscore >= thresh, 
         Network == "STRING") %>% 
  mutate(prop_new_genes =
           nNewGenes/nSeeds)

string %>%
  ggplot(aes(x = Features, 
             y = prop_new_genes)) +
  geom_jitter(shape=16, 
              colour = "#90A4AE",
              position=position_jitter(0.2)) +
  geom_boxplot(color = c("#C62828"),
               fill = "NA",
               outlier.size=0, 
               outlier.shape=NA) +
  theme_classic()

ggsave(file = paste0(out_dir, "/string_proportion_new_genes_boxplots.pdf"),
       width = 3,
       height = 3)


# find and move the diseases with cv >= thresh to new folder --------------
# find diseases with SLA cv score >= thresh
# and save them to refer to later
load("../results/GenePlexus_parameters/Geneplexus_summary.Rdata")

sla1 = 
  summary_df %>%
  filter(Network == "STRING",
         Features == "Adjacency",
         CVscore >= thresh) %>%
  pull(Disease)

remove_trait = c("chronic_inflammation_gene_shot",
                 "chronic_inflammation_gene_shot_pubs_greater10",
                 "chronic_inflammation_go",
                 "Height_gene_shot")

sla1 = sla1[!sla1 %in% remove_trait]
save(sla1, file = "../results/GenePlexus_parameters/diseases_with_SLA_models_cv_greater1.Rdata")

# cp predictions with cv > thresh into new folder
print(summary_df$CVscore)
to_copy = 
  summary_df %>%
  filter(Network == "STRING",
         Features == "Adjacency",
         CVscore >= thresh,
         !Disease %in% c("chronic_inflammation_gene_shot",
                         "chronic_inflammation_gene_shot_pubs_greater10",
                         "Height_gene_shot")) %>%
  mutate(to_copy = 
           paste0(
             Disease,
             "--",
             Network,
             "--",
             Features,
             "--",
             negativesFrom, 
             "--predictions.tsv")) %>%
  pull(to_copy)

new_dir = paste0("../GenePlexus_output/predictions_cv_greater", thresh)
if(!dir.exists(new_dir)){dir.create(new_dir)}

file.copy(paste0(args[1], to_copy), new_dir)
warnings()

