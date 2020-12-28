library(tidyverse)
library(ggsci)
library(cowplot)
library(gridExtra)
library(grid)
library(ggpubr)

get.tables.dir.path <- function(){
  return("tables")
}

clin <- read_csv(file.path(get.tables.dir.path(), "clinical_multi_omics.csv"))
surv <- read_csv(file.path(get.tables.dir.path(), "survival_multi_omics.csv"))

all_surv <- surv %>% mutate(tool=X1) %>% select(-X1) %>% gather(cancer, survival, -tool)
all_clin <- clin %>% mutate(tool=X1) %>% select(-X1) %>% gather(cancer, clin, -tool)

all_clin$is_sumo <- grepl("SUMO",all_clin$tool)
sumo_clin <- all_clin %>% filter(is_sumo == TRUE) %>% separate(tool, c("tmp", "k"), sep = "SUMO") %>% select(-tmp) %>% mutate(tool=paste0("SUMO",k))

all_surv$is_sumo <- grepl("SUMO",all_surv$tool)
sumo_surv <- all_surv %>% filter(is_sumo == TRUE) %>% separate(tool, c("tmp", "k"), sep = "SUMO") %>% select(-tmp) %>% mutate(tool=paste0("SUMO",k))

sumo_data <- sumo_surv %>% full_join(sumo_clin) %>% filter(k != "")
sumo_data

cancers <- c("AML", "BIC", "COAD", "GBM", "KIRC", "LIHC", "LUSC", "SKCM", "OV", "SARC")
selected <- c(12, 13, 7, 10, 6, 2, 9, 7, 9, 6)
selected_sumo <- sumo_data %>% inner_join(tibble(k=as.character(selected), cancer=cancers)) %>% mutate(tool="SUMO")
selected_sumo

clin_data <- all_clin %>% filter(is_sumo == FALSE)
surv_data <- all_surv %>% filter(is_sumo == FALSE)

data <- clin_data %>% full_join(surv_data) %>% full_join(selected_sumo %>% select(-k))
data %>% filter(tool == "iClusterBayes")

###
data <- read_tsv("benchmark_results_selected_0.2.5.tsv")
data <- data %>% mutate(tool = as.factor(tool), is_sumo = as.logical(is_sumo))
data$tool <- factor(data$tool, levels = c("NEMO", "SUMO","LRAcluster", "MCCA", 
                                          "PINSPlus", "SNF", "iClusterBayes"))

data %>% mutate(sig_surv = survival >= -log10(0.05), sig_clin = clin > 0) %>% filter(sig_surv == TRUE | sig_clin == TRUE) %>%
  group_by(tool) %>% summarise(clinical=sum(sig_clin), survival=sum(sig_surv)) %>%
  select(tool, survival, clinical)
stats <- .Last.value
stats <- tibble(tool=c('LRAcluster', 'MCCA', 'NEMO', 'PINSPlus', 'SNF', 'iClusterBayes', 'SUMO')) %>%
  mutate(tool=as.factor(tool)) %>%
  left_join(stats)

colnames(stats) <- c("Method", "Number of cancers\n with differential survival", "Number of cancers\n with clinical enrichment")

tt <- ttheme_default(core=list(
  fg_params=list(fontface=c(rep("plain", 6), "bold")),
  bg_params = list(fill=c("grey95", "grey90","grey95", "grey90","grey95", "grey90","grey60"),
                   alpha = rep(c(1,0.5), each=6))
), base_size = 10.5)

figB <- tableGrob(stats, rows=NULL, theme=tt)

ggplot(data, aes(survival, clin, color=tool, shape=is_sumo)) + 
  geom_point(size=3) + 
  scale_y_continuous(breaks=c(0,1,2,3,4), 
                     labels=c("0","1","2","3","4"),
                     limits=c(0,NA)) +
  xlab(expression(-log[10](`logrank p-value`))) + 
  ylab("number of enriched clinical parameters") +
  facet_wrap(cancer ~ ., ncol=4) +
  geom_vline(xintercept=-log10(0.05), color="red") + 
  theme_bw() + 
  theme(legend.position=c(0.75,0.15), 
        legend.direction = "horizontal",
        panel.grid.minor.y = element_blank()) +
  guides(shape=FALSE, 
         color=guide_legend(override.aes = list(shape=c(16,17,16,16,16,16,16), size=3), nrow=3, byrow = TRUE)) +
  labs(color="", shape="") +
  scale_color_npg() 

figA <- last_plot()

cairo_pdf("benchmark_results.pdf", width=7, height=7.5)
ggarrange(figA, figB, labels = c('A','B'), ncol=1, heights = c(2,1))
dev.off()
