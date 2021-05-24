library(tidyverse)
library(ggsci)
library(mclust)
library(ggpubr)

get.tables.dir.path <-function(){
  return(file.path("results_eval", paste0('results_', REP,'/tables')))
}

get.tool.labels <- function(tool, rep){
  return(file.path("results_eval", paste0('results_', rep), paste("gbm",tool, "all", sep = "_")))
}

clin_all <- NULL
surv_all <- NULL
for (REP in 1:10){
  clin <- read_csv(file.path(get.tables.dir.path(), "clinical_multi_omics.csv"), col_types = cols())
  surv <- read_csv(file.path(get.tables.dir.path(), "survival_multi_omics.csv"), col_types = cols())
  colnames(clin) <- c('X1', paste0("GBM_", REP))
  colnames(surv) <- c('X1', paste0("GBM_", REP))
  if (REP == 1){
    clin_all <- clin
    surv_all <- surv
  } else {
    clin_all <- full_join(clin_all, clin)
    surv_all <- full_join(surv_all, surv)
  }
}

all_surv <- surv_all %>% mutate(tool=X1) %>% select(-X1) %>% gather(cancer, survival, -tool)
all_clin <- clin_all %>% mutate(tool=X1) %>% select(-X1) %>% gather(cancer, clin, -tool)
data <- all_surv %>% full_join(all_clin, by = c("tool", "cancer"))

SUMO_k <- 12
data_other <- data %>% filter(!grepl("SUMO", tool))
data <- data %>% 
  filter(tool == paste0("SUMO", SUMO_k)) %>%
  mutate(tool = "SUMO") %>%
  full_join(data_other,  by = c("tool", "cancer", "survival", "clin"))

data

data[data$tool == "PINS",]$tool = "PINSPlus"

tools <- c("nemo", paste0("sumo", SUMO_k), "lracluster", "mcca", "pins", "snf", "iCluster", "cimlr")
names(tools) <- c("NEMO", "SUMO","LRAcluster", "MCCA", "PINSPlus", "SNF", "iClusterBayes", "CIMLR")

data <- data %>% mutate(tool = as.factor(tool))
data$tool <- factor(data$tool, levels = names(tools))

tool_pal <- pal_npg("nrc")(10)[c(1:7,10)]
names(tool_pal) <- names(tools)

data <- data %>% separate(cancer, c("cancer_type", "idx"), sep="_", remove = FALSE)
data$tool_name <- sapply(data$tool, function(x){tools[names(tools) == x]})
data <- data %>% mutate(path = get.tool.labels(tool_name, idx))

ARI_data <- NULL
for (the_tool in data$tool %>% unique()){
  paths <- data %>% filter(tool==the_tool) %>% pull(path)
  labels <- NULL
  for (fname in paths){
    load(fname)
    results_rep <- strsplit(fname, "/")[[1]][2]
    if (is.null(labels)){
      labels <- tibble(!!(results_rep) := clustering)
    } else {
      labels <- add_column(labels, !!(results_rep) := clustering)
    }
  }
  
  ari <- sapply(1:ncol(labels), function(x){vec1 <- labels %>% pull(!!(x)); 
                                            sapply(1:ncol(labels), function(y){adjustedRandIndex(vec1, labels %>% pull(!!(y)))})})
  ARI_tbl <- tibble(tool=the_tool, medianARI=median(ari[upper.tri(ari, diag=FALSE)]), sdARI=sd(ari[upper.tri(ari, diag=FALSE)]))
  if (is.null(ARI_data)){
    ARI_data <- ARI_tbl
  } else {
    ARI_data <- ARI_data %>% full_join(ARI_tbl, by = c("tool", "medianARI", "sdARI"))
  }
}

ARI_data <- ARI_data %>% mutate(tool= as.factor(tool)) 
ARI_data$tool <- factor(ARI_data$tool, levels = names(tools))

ARI_data$tool %>% unique()
ARI_data %>%
  ggplot() + geom_point(aes(x=tool, color=tool, y=medianARI), size=2) +
  geom_errorbar(aes(x = tool, ymin=medianARI-sdARI, ymax=medianARI+sdARI, color=tool), size=1, width=0.4) + 
  geom_text(aes(x=tool, y=c(0.8, 0.9, 0.87, 0.88, 0.4, 0.85, 0.74, 0.88), label=tool)) +
  scale_color_manual(values=tool_pal) +
  theme_bw() + 
  theme(legend.position = "null", axis.text = element_text(size=12),
        strip.text.x = element_text(size = 12, face="bold"), axis.title = element_text(size=12), 
        legend.text=element_text(size=12), axis.text.y = element_blank(), axis.ticks.y = element_blank()) +
  labs(x="",y="median pairwise ARI") +
  coord_flip()
stability <- last_plot()


pdf("evaluate_benchmark.pdf", width = 11, height = 7)

p <-  data %>%
    ggplot() + geom_point(aes(x=survival,y=clin,color=as.factor(tool)), size=3) + facet_wrap(tool~., nrow = 3) +
    theme_bw() + labs(color="tool") + geom_vline(xintercept = -log10(0.05), color="red") + 
    theme(legend.position = "null", axis.text = element_text(size=12), strip.text.x = element_text(size = 12, face="bold"), 
          axis.title = element_text(size=12), legend.text=element_text(size=12)) + scale_color_manual(values=tool_pal) +
    ylim(0,2) + xlim(0.9, 5) + ylab('# enriched clinical parameters') + xlab("-log10(logrank pvalue)") 

ggarrange(p, stability, ncol=2, widths = c(3.5, 1.5),labels = c('A','B'))

dev.off()
