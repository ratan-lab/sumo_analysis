library(tidyverse)
library(Rtsne)
library(ggpubr)
setwd("~/Desktop/sumo_analysis/simulations/noisy")
  
map_vals <- function(val, sampling = "gamma"){
  stopifnot(sampling %in% c('gamma', 'gauss'))
  if (sampling == 'gamma'){
    RANGE <- gamma_spl_k
  } else {
    RANGE <- gauss_spl_std
  }
  return(RANGE[as.integer(val)])
}

gamma <- read_tsv("gamma.tsv")
gamma <- gamma %>% separate(subtype, c("sampling", "var")) %>% mutate(var = map_vals(var, sampling = 'gamma'))
gamma %>% ggplot() + geom_line(aes(x=var,y=val, group=algoritm, color=algoritm)) + 
  facet_wrap(metric~.) + theme(legend.position = "bottom") + ylab("metric value") + xlab("k") + 
  ggtitle(paste0('Gamma(theta=',gamma_spl_theta,') & Gauss(mean=',gamma_spl_mean,', sd=',gamma_spl_std,')'))

gamma %>% filter(metric=='ARI') %>% ggplot() + geom_line(aes(x=var,y=val, group=algoritm, color=algoritm)) + 
  facet_wrap(algoritm~., ncol=4) + theme(legend.position = "bottom") + ylab("metric value") + xlab("k") + 
  ggtitle(paste0('Gamma(theta=',gamma_spl_theta,') & Gauss(mean=',gamma_spl_mean,', sd=',gamma_spl_std,')'))


gauss <- read_tsv("gauss.tsv")
gauss <- gauss %>% separate(subtype, c("sampling", "var")) %>% mutate(var = map_vals(var, sampling = 'gauss'))
gauss %>% ggplot() + geom_line(aes(x=var,y=val, group=algoritm, color=algoritm)) + 
  facet_wrap(metric~.) + theme(legend.position = "bottom") + ylab("metric value") + xlab("mean") +
  ggtitle(paste0('Gauss(mean=',gauss_spl_mean,') & Gamma(k=',gauss_spl_k,', theta=',gauss_spl_theta,')'))

gauss %>% filter(metric=='ARI') %>% ggplot() + geom_line(aes(x=var,y=val, group=algoritm, color=algoritm)) + 
  facet_wrap(algoritm~., ncol=4) + theme(legend.position = "bottom") + ylab("metric value") + xlab("mean") +
  ggtitle(paste0('Gauss(mean=',gauss_spl_mean,') & Gamma(k=',gauss_spl_k,', theta=',gauss_spl_theta,')'))
