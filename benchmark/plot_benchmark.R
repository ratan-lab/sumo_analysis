library(tidyverse)

clin <- read_csv(file.path(get.tables.dir.path(), "clinical_multi_omics.csv"))
surv <- read_csv(file.path(get.tables.dir.path(), "survival_multi_omics.csv"))

###

all_surv <- read_csv(surv_fname) %>% mutate(tool=X1) %>% select(-X1) %>% gather(cancer, survival, -tool)
all_clin <- read_csv(clin_fname) %>% mutate(tool=X1) %>% select(-X1) %>% gather(cancer, clin, -tool)

all_clin$is_sumo <- grepl("SUMO",all_clin$tool)
sumo_clin <- all_clin %>% filter(is_sumo == TRUE) %>% separate(tool, c("tmp", "k"), sep = "SUMO") %>% select(-tmp) %>% mutate(tool=paste0("SUMO",k))

all_surv$is_sumo <- grepl("SUMO",all_surv$tool)
sumo_surv <- all_surv %>% filter(is_sumo == TRUE) %>% separate(tool, c("tmp", "k"), sep = "SUMO") %>% select(-tmp) %>% mutate(tool=paste0("SUMO",k))

sumo_data <- sumo_surv %>% full_join(sumo_clin) %>% filter(k != "")

sumo_data

sumo_data %>%
  mutate(tool="SUMO") %>%
  ggplot() + geom_point(aes(x=survival,y=clin,color=as.factor(k)), size=2) + facet_wrap(cancer~., scales = "free") + 
  theme(legend.position = "bottom") + ggtitle("SUMO") + labs(color="k") + geom_vline(xintercept = -log10(0.05), color="red")

###

clin_long <- clin %>%
  gather(cancer, clin, -X1) %>%
  mutate(tool=X1) %>%
  select(-X1)

surv %>%
  gather(cancer, surv, -X1) %>%
  mutate(tool=X1) %>%
  select(-X1) %>%
  full_join(clin_long) %>%
  ggplot() + geom_point(aes(x=surv, y=clin, color=tool)) + facet_wrap(cancer~., scales="free") +
  xlab("-log10(logrank pvalue)") + ylab("# enriched clinical parameters") + labs(color="") +
  theme(legend.position="bottom") + geom_vline(xintercept = -log10(0.05), color="red")

p <- last_plot()