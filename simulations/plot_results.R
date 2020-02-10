library(tidyverse)
library(ggsci)

gauss_spl_std <- seq(0,4, 0.2)

map_vals <- function(val){
  RANGE <- gauss_spl_std
  return(RANGE[as.integer(val)])
}

data <- read_tsv("double_gauss.tsv")
data <- data %>% separate(subtype, c("a", "b", "var")) %>% mutate(var = map_vals(var)) %>% select(-a, -b)
dataARI <- data %>% filter(metric == "ARI")

dataARI <- dataARI %>% group_by(algoritm, var) %>% 
  summarise(medianARI= median(val), minARI=min(val), maxARI=max(val)) %>% full_join(dataARI)

dataARI %>%
  ggplot() + 
    geom_ribbon(aes(x = var, ymin = minARI, ymax = maxARI, fill=algoritm), alpha=0.3) +
    geom_line(aes(x=var,y=medianARI, group=algoritm, color=algoritm), size=1) + facet_wrap(algoritm~.) +
    theme_bw() + ylab("median ARI") + xlab("σ") + 
    theme(legend.position = "bottom", axis.text = element_text(size=12), strip.text.x = element_text(size = 12, face="bold"), 
          axis.title = element_text(size=12), legend.text=element_text(size=12), legend.title = element_blank()) +
    scale_color_npg() + scale_fill_npg()

dataARI %>%
  ggplot() + 
  geom_line(aes(x=var,y=medianARI, group=algoritm, color=algoritm), size=1) + theme_bw() + ylab("median ARI") + xlab("σ") + 
  theme(legend.position = "bottom", axis.text = element_text(size=12), strip.text.x = element_text(size = 12, face="bold"), 
        axis.title = element_text(size=12), legend.text=element_text(size=12), legend.title = element_blank()) +
  scale_color_npg()
