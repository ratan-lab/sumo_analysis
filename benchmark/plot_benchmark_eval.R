library(tidyverse)
library(ggsci)

get.tables.dir.path <-function(){
  return(paste0('results_', REP,'/tables'))
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
data <- all_surv %>% full_join(all_clin)
data

data %>%
  ggplot() + geom_point(aes(x=survival,y=clin,color=as.factor(tool)), size=3) + facet_wrap(cancer~., ncol=4) + 
  theme_bw() + labs(color="tool") + geom_vline(xintercept = -log10(0.05), color="red") + ylab('# enriched clinical parameters') +
  xlab("-log10(logrank pvalue)") + theme(legend.position = "bottom", axis.text = element_text(size=12),
                                         strip.text.x = element_text(size = 12, face="bold"), axis.title = element_text(size=12), 
                                         legend.text=element_text(size=12)) +scale_color_npg()

data %>%
  ggplot() + geom_point(aes(x=survival,y=clin,color=as.factor(tool)), size=3) + facet_wrap(tool~.) + 
  theme_bw() + labs(color="tool") + geom_vline(xintercept = -log10(0.05), color="red") + ylab('# enriched clinical parameters') +
  xlab("-log10(logrank pvalue)") + theme(legend.position = "null", axis.text = element_text(size=12),
                                         strip.text.x = element_text(size = 12, face="bold"), axis.title = element_text(size=12), 
                                         legend.text=element_text(size=12)) +scale_color_npg()
