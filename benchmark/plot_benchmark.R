library(tidyverse)
library(ggsci)

clin <- read_csv(file.path(get.tables.dir.path(), "clinical_multi_omics.csv"))
surv <- read_csv(file.path(get.tables.dir.path(), "survival_multi_omics.csv"))

###

all_surv <- surv %>% mutate(tool=X1) %>% select(-X1) %>% gather(cancer, survival, -tool)
all_clin <- clin %>% mutate(tool=X1) %>% select(-X1) %>% gather(cancer, clin, -tool)

all_clin$is_sumo <- grepl("SUMO",all_clin$tool)
sumo_clin <- all_clin %>% filter(is_sumo == TRUE) %>% separate(tool, c("tmp", "k"), sep = "SUMO") %>% select(-tmp) %>% mutate(tool=paste0("SUMO",k))

all_surv$is_sumo <- grepl("SUMO",all_surv$tool)
sumo_surv <- all_surv %>% filter(is_sumo == TRUE) %>% separate(tool, c("tmp", "k"), sep = "SUMO") %>% select(-tmp) %>% mutate(tool=paste0("SUMO",k))

sumo_data <- sumo_surv %>% full_join(sumo_clin) %>% filter(k != "")

sumo_data

sumo_data %>%
  mutate(tool="SUMO") %>%
  ggplot() + geom_point(aes(x=survival,y=clin,color=as.factor(k)), size=3) + facet_wrap(cancer~., scales = "free", ncol=4) + 
  theme_bw() + labs(color="k") + geom_vline(xintercept = -log10(0.05), color="red") + ylab('# enriched clinical parameters') +
  xlab("-log10(logrank pvalue)") + theme(legend.position = "bottom", axis.text = element_text(size=12),
        strip.text.x = element_text(size = 12, face="bold"), axis.title = element_text(size=12), 
        legend.text=element_text(size=12))

sumo_plot <- last_plot()

###
cancers <- c("AML", "BIC", "COAD", "GBM", "KIRC", "LIHC", "LUSC", "SKCM", "OV", "SARC")
selected <- c(12, 13, 7, 10, 6, 2, 9, 7, 9, 6)

selected_sumo <- sumo_data %>% inner_join(tibble(k=as.character(selected), cancer=cancers)) %>% mutate(tool="SUMO")
selected_sumo

clin_data <- all_clin %>% filter(is_sumo == FALSE)
surv_data <- all_surv %>% filter(is_sumo == FALSE)

data <- clin_data %>% full_join(surv_data) %>% full_join(selected_sumo %>% select(-k))
data
data %>%
  ggplot() + geom_point(aes(x=survival, y=clin, color=tool, shape= is_sumo), size = 3) + 
  facet_wrap(cancer~., scales="free", ncol=4) + xlab("-log10(logrank pvalue)") + ylab("# enriched clinical parameters") + 
  labs(color="") + theme_bw() + geom_vline(xintercept = -log10(0.05), color="red") + scale_shape(guide="none") +
  theme(legend.position = "bottom", axis.text = element_text(size=12), strip.text.x = element_text(size = 12, face="bold"), 
        axis.title = element_text(size=12), legend.text=element_text(size=12)) + scale_color_npg()

results_plot <- last_plot()

###

data %>% mutate(sig_surv = survival >= -log10(0.05), sig_clin = clin > 0) %>% filter(sig_surv == TRUE | sig_clin == TRUE) %>%
  group_by(tool) %>% summarise(clinincal=sum(sig_clin), survival=sum(sig_surv))
