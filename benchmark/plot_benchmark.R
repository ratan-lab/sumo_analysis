library(tidyverse)
library(ggsci)
library(cowplot)
library(gridExtra)
library(grid)
library(ggpubr)

get.tables.dir.path <- function(){
  return("results/tables")
}

clin <- read_csv(file.path(get.tables.dir.path(), "clinical_multi_omics.csv"))
surv <- read_csv(file.path(get.tables.dir.path(), "survival_multi_omics.csv"))
pathways <- read_tsv("benchmark_pathway_activity.tsv") %>% select(-total)
pathways %>% 
  gather(pathway, pval, -tool,-cancer) %>%
  filter(pval <= 0.05) %>%
  group_by(tool, cancer) %>%
  summarise(n=n()) %>%
  filter(n >0) %>%
  group_by(tool) %>%
  summarise(n=n())
pathway_stats <- .Last.value

all_surv <- surv %>% mutate(tool=X1) %>% select(-X1) %>% gather(cancer, survival, -tool)
all_clin <- clin %>% mutate(tool=X1) %>% select(-X1) %>% gather(cancer, clin, -tool)

all_clin$is_sumo <- grepl("SUMO",all_clin$tool)
sumo_clin <- all_clin %>% filter(is_sumo == TRUE) %>% separate(tool, c("tmp", "k"), sep = "SUMO") %>% select(-tmp) %>% mutate(tool=paste0("SUMO",k))

all_surv$is_sumo <- grepl("SUMO",all_surv$tool)
sumo_surv <- all_surv %>% filter(is_sumo == TRUE) %>% separate(tool, c("tmp", "k"), sep = "SUMO") %>% select(-tmp) %>% mutate(tool=paste0("SUMO",k))

sumo_data <- sumo_surv %>% full_join(sumo_clin) %>% filter(k != "")
sumo_data

cancers <- c("AML", "BIC", "COAD", "GBM", "KIRC", "LIHC", "LUSC", "SKCM", "OV", "SARC")
selected <- c(11, 2, 2, 12, 6, 10, 12, 11, 3, 6) 

selected_sumo <- sumo_data %>% inner_join(tibble(k=as.character(selected), cancer=cancers)) %>% mutate(tool="SUMO")
selected_sumo

clin_data <- all_clin %>% filter(is_sumo == FALSE)
surv_data <- all_surv %>% filter(is_sumo == FALSE)

data <- clin_data %>% full_join(surv_data) %>% full_join(selected_sumo %>% select(-k))

data[data$tool == "PINS",]$tool = "PINSPlus"
data <- data %>% mutate(tool = as.factor(tool), is_sumo = as.logical(is_sumo))
data$tool <- factor(data$tool, levels = c("NEMO", "SUMO","LRAcluster", "MCCA", 
                                          "PINSPlus", "SNF", "iClusterBayes", "CIMLR"))

data %>% mutate(sig_surv = survival >= -log10(0.05), sig_clin = clin > 0) %>% filter(sig_surv == TRUE | sig_clin == TRUE) %>%
  group_by(tool) %>% summarise(clinical=sum(sig_clin), survival=sum(sig_surv)) %>%
  select(tool, survival, clinical) %>% left_join(pathway_stats) 
stats <- .Last.value
stats <- tibble(tool=c('LRAcluster', 'MCCA', 'NEMO', 'PINSPlus', 'SNF', 'iClusterBayes', 'CIMLR', 'SUMO')) %>%
  mutate(tool=as.factor(tool)) %>%
  left_join(stats)

colnames(stats) <- c("Method", "Number of cancers\nwith differential survival", "Number of cancers\nwith clinical enrichment",
                     "Number of cancers\nwith pathway activity\nenrichment")

tt <- ttheme_default(core=list(
  fg_params=list(fontface=c(rep("plain", 7), "bold")),
  bg_params = list(fill=c("grey95", "grey90","grey95", "grey90","grey95", "grey90","grey95","grey60"),
                   alpha = rep(c(1,0.5), each=6))
), base_size = 10.5)

figB <- tableGrob(stats, rows=NULL, theme=tt)

tool_pal <- pal_npg("nrc")(10)[c(1:7,10)]
names(tool_pal) <- c("NEMO", "SUMO","LRAcluster", "MCCA", "PINSPlus", "SNF", "iClusterBayes", "CIMLR")

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
         color=guide_legend(override.aes = list(shape=c(16,17,16,16,16,16,16,16), size=3), nrow=3, byrow = TRUE)) +
  labs(color="", shape="") +
  scale_color_manual(values=tool_pal)
  
figA <- last_plot()

pdf("benchmark_results.pdf", width=7, height=9)
ggarrange(figA, figB, labels = c('A','B'), ncol=1, heights = c(2,1))
dev.off()
  