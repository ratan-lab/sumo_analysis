library(tidyverse)
setwd("~/Desktop/sumo_analysis/simulations/noisy")

map_vals <- function(val){
  RANGE <- gauss_spl_std
  return(RANGE[as.integer(val)])
}

data <- read_tsv("double_gauss.tsv")
data <- data %>% separate(subtype, c("sampling", "tmp", "var")) %>% mutate(var = map_vals(var), sampling= paste(sampling, tmp, sep="_")) %>% select(-tmp)
data %>% ggplot() + geom_line(aes(x=var,y=val, group=algoritm, color=algoritm)) + 
  facet_wrap(metric~.) + theme(legend.position = "bottom") + ylab("metric value") + xlab("k") + 
  ggtitle(paste0('Gauss(mean=',gauss_spl_mean,') & Gauss(mean=',gauss_mean,', sd=',gauss_std,')'))

data %>% filter(metric == 'ARI') %>% ggplot() + geom_line(aes(x=var,y=val, group=algoritm, color=algoritm)) + 
  facet_wrap(algoritm~., ncol=4) + theme(legend.position = "null") + ylab("ARI") + xlab("k") + 
  ggtitle(paste0('Gauss(mean=',gauss_spl_mean,') & Gauss(mean=',gauss_mean,', sd=',gauss_std,')'))

