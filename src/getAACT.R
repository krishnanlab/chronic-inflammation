# download clinicaltrials.gov data
#-------------
# mkdir ../data/clinical_trials
# mkdir ../data/clinical_trials/20211004
# cd ../data/clinical_trials/20211004
# wget https://aact.ctti-clinicaltrials.org/static/exported_files/monthly/20211004_pipe-delimited-export.zip
#-------------
# use instructions from https://aact.ctti-clinicaltrials.org/pipe_files_with_r

# summary info
studies = read.table(file =  "studies.txt",
                     header =  TRUE,
                     sep =  "|",
                     na.strings =  "",
                     stringsAsFactors =  FALSE, 
                     comment.char =  "",
                     quote =  "\"",
                     fill =  FALSE)

studies = studies %>%
  select(nct_id,
         phase,
         overall_status,
         completion_date)

#diseases
# mesh terms for diseases
browse_conditions = read.table(file =  "browse_conditions.txt",
                               header =  TRUE,
                               sep =  "|",
                               na.strings =  "",
                               stringsAsFactors =  FALSE, 
                               comment.char =  "",
                               quote =  "\"",
                               fill =  FALSE)

colnames(browse_conditions)[3] = "condition_mesh_term"
browse_conditions$id = NULL
browse_conditions$downcase_mesh_term = NULL

# drugs live here
# only ones with mesh ids
browse_interventions = read.table(file =  "browse_interventions.txt",
                                  header =  TRUE,
                                  sep =  "|",
                                  na.strings =  "",
                                  stringsAsFactors =  FALSE, 
                                  comment.char =  "",
                                  quote =  "\"",
                                  fill =  FALSE)

colnames(browse_interventions)[3] = "intervention_mesh_term"
browse_interventions$id = NULL
browse_interventions$downcase_mesh_term = NULL

# bind together
cond_int = left_join(browse_conditions, 
                     browse_interventions, 
                     by = "nct_id")

ct = left_join(cond_int, studies, by = "nct_id")
save(ct, file = "../data/clinical_trials/20211004/basic_clinical_trial_info.rdata")