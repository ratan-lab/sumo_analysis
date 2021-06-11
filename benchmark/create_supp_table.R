# creates Supplementary Table S2
library(tidyverse)

# OTHER TOOLS
# check survival
surv_table <- read_csv("results/tables/survival_multi_omics.csv") %>% mutate(tool = X1) %>% select(-X1) %>% gather(cancer, surv_pval, -tool)
surv_table[surv_table$tool == "PINS",]$tool <- "PINSPlus"
surv_table_other <- surv_table %>% filter(!grepl("SUMO", tool))
# check clinical
clin_table <- read_csv("results/tables/clinical_multi_omics.csv") %>% mutate(tool = X1) %>% select(-X1) %>% gather(cancer, clin_params, -tool)
clin_table[clin_table$tool == "PINS",]$tool <- "PINSPlus"
clin_table_other <- clin_table %>% filter(!grepl("SUMO", tool))
# check runtime
time_table <- read_csv("results/tables/runtime_multi_omics.csv") %>% mutate(tool = X1) %>% select(-X1) %>% gather(cancer, runtime, -tool)
time_table[time_table$tool == "PINS",]$tool <- "PINSPlus"
time_table_other <- time_table %>% filter(!grepl("SUMO", tool))

# SUMO 
# check survival, clinical & runtime
surv_table_sumo <- surv_table %>% filter(grepl("SUMO", tool)) %>% separate(tool, c('tmp', 'k'), sep="SUMO") %>% 
  select(-tmp) %>% filter(k != "") %>% mutate(tool = "SUMO", k= as.numeric(k))
clin_table_sumo <- clin_table %>% filter(grepl("SUMO", tool)) %>% separate(tool, c('tmp', 'k'), sep="SUMO") %>% 
  select(-tmp) %>% filter(k != "") %>% mutate(tool = "SUMO", k= as.numeric(k))
time_table_sumo <- time_table %>% filter(grepl("SUMO", tool)) %>% separate(tool, c('tmp', 'k'), sep="SUMO") %>% 
  select(-tmp) %>% filter(k != "") %>% mutate(tool = "SUMO", k= as.numeric(k))

cancers <- c("AML", "BIC", "COAD", "GBM", "KIRC", "LIHC", "LUSC", "SKCM", "OV", "SARC")
selected <- c(12, 13, 7, 10, 6, 2, 9, 7, 9, 6)
surv_table_sumo <- tibble(cancer=cancers, k=selected) %>% left_join(surv_table_sumo, by = c("cancer", "k")) %>% 
  select(tool, cancer, surv_pval)
clin_table_sumo <- tibble(cancer=cancers, k=selected) %>% left_join(clin_table_sumo, by = c("cancer", "k")) %>% 
  select(tool, cancer, clin_params)
time_table_sumo <- tibble(cancer=cancers, k=selected) %>% left_join(time_table_sumo, by = c("cancer", "k")) %>% 
  select(tool, cancer, runtime)

###
surv_table_final <- surv_table_sumo %>% full_join(surv_table_other, by = c("tool", "cancer", "surv_pval"))
clin_table_final <- clin_table_sumo %>% full_join(clin_table_other, by = c("tool", "cancer", "clin_params"))
time_table_final <- time_table_sumo %>% full_join(time_table_other, by = c("tool", "cancer", "runtime"))

#num clusters
nclus_table <- read_csv("results/tables/num_cluster_multi_omics.csv") %>% mutate(tool = X1) %>% select(-X1) %>% gather(cancer, k, -tool)
nclus_table[nclus_table$tool == "PINS",]$tool <- "PINSPlus"
nclus_table_final <- nclus_table %>% filter(!grepl("SUMO", tool)) %>% 
  full_join(tibble(tool="SUMO", cancer=cancers, k=selected), by = c("tool", "cancer", "k"))

res_table <- surv_table_final %>% left_join(clin_table_final, by = c("tool", "cancer")) %>% 
  left_join(nclus_table_final, by = c("tool", "cancer"))

# which clinical params
fnames <- dir('results')[grepl('*_clin', dir('results'))]
cancer_codes <-  c("aml", "breast", "colon", "gbm", "kidney", "liver", "lung", "melanoma", "ovarian", "sarcoma")
names(cancer_codes) <- c('AML', 'BIC', 'COAD', 'GBM', 'KIRC', 'LIHC', 'LUSC', 'SKCM', 'OV', 'SARC')
clins <- NULL
for (fname in fnames){
  load(file.path('results',fname))
  cancer_name <- strsplit(fname, "_")[[1]][1]
  tool_name <- strsplit(fname, "_")[[1]][2]
  tib <- tibble(cancer=names(cancer_codes[cancer_codes==cancer_name]), tool=tool_name, param=names(enrichment.pvalues), pval=enrichment.pvalues)
  if (is.null(clins)){
    clins <- tib
  } else {
    clins <- clins %>% full_join(tib, by = c("cancer", "tool", "param", "pval"))
  }
}

tool_names <- c("iCluster", "lracluster", "mcca", "nemo", "pins", "snf", "cimlr") 
names(tool_names) <- c("iClusterBayes", "LRAcluster", "MCCA", "NEMO", "PINSPlus", "SNF", "CIMLR")
clins_others <- clins %>% filter(!grepl("sumo", tool))
clins_others$tool <- unlist(lapply(clins_others$tool, function(x){names(tool_names[tool_names == x])}))
clins_others <- clins_others %>% spread(param, pval)

clins_sumo <- clins %>% filter(grepl("sumo", tool))
clins_sumo <- clins_sumo %>% separate(tool, c('tmp', 'k'), sep="sumo") %>% select(-tmp) %>% filter(k != "") %>% 
  mutate(tool = "SUMO", k=as.numeric(k)) %>% spread(param, pval)
clins_sumo <- tibble(tool="SUMO", cancer=cancers, k=selected) %>% left_join(clins_sumo, by = c("tool", "cancer", "k"))

clins_full <- clins_others %>% full_join(clins_sumo  %>% select(-k),
                                         by = c("cancer", "tool", "age_at_initial_pathologic_diagnosis", "gender", 
                                                "pathologic_M", "pathologic_N", "pathologic_stage", "pathologic_T"))
clins_full
clins_full  %>% 
  mutate(a = !is.na(age_at_initial_pathologic_diagnosis), b = !is.na(gender), c = !is.na(pathologic_M),
         d = !is.na(pathologic_N), e = !is.na(pathologic_stage), f = !is.na(pathologic_T)) %>%
  group_by(tool, cancer) %>%
  mutate(tested_params = sum(a,b,c,d,e,f)) %>%
  select(-a,-b,-c,-d,-e,-f) %>%
  mutate(age_at_initial_pathologic_diagnosis = age_at_initial_pathologic_diagnosis * tested_params,
         gender = gender * tested_params, pathologic_M = pathologic_M * tested_params,
         pathologic_N = pathologic_N * tested_params, pathologic_stage = pathologic_stage * tested_params,
         pathologic_T = pathologic_T * tested_params) %>%
  mutate(age_at_initial_pathologic_diagnosis = ifelse(age_at_initial_pathologic_diagnosis < 0.05, 1, 0),
         gender = ifelse(gender < 0.05, 1, 0),
         pathologic_M = ifelse(pathologic_M < 0.05, 1, 0),
         pathologic_N = ifelse(pathologic_N < 0.05, 1, 0),
         pathologic_stage = ifelse(pathologic_stage < 0.05, 1, 0),
         pathologic_T = ifelse(pathologic_T < 0.05, 1, 0)) %>%
  mutate(sum_clin_params=sum(age_at_initial_pathologic_diagnosis, gender, pathologic_M, pathologic_N, pathologic_stage, pathologic_T, na.rm = TRUE))

clins_summary <- .Last.value
clins_summary %>% select(cancer, tool, sum_clin_params)  %>% left_join(res_table, by = c("cancer", "tool")) %>% 
  filter(sum_clin_params != clin_params) #-> no diff

res_table <- res_table %>% left_join(clins_summary %>% select(-sum_clin_params, -tested_params), by = c("tool", "cancer"))
res_table <- res_table %>% left_join(time_table_final, by=c('tool', 'cancer'))

write_tsv(res_table %>% arrange(cancer, tool), "benchmark_results_selected_0.2.6.tsv")