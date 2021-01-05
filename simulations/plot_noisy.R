library(tidyverse)
library(ggsci)

gauss_spl_std <- seq(0,4, 0.2)

map_vals <- function(val){
  RANGE <- gauss_spl_std
  return(RANGE[as.integer(val)])
}

data <- read_tsv("noisy/double_gauss.tsv")
data <- data %>% separate(subtype, c("a", "b", "var")) %>% mutate(var = map_vals(var)) %>% select(-a, -b)
dataARI <- data %>% filter(metric == "ARI")

dataARI <- dataARI %>% group_by(tool, var) %>% 
  summarise(medianARI= median(val), minARI=min(val), maxARI=max(val)) %>% full_join(dataARI)

dataARI <- dataARI %>% ungroup() %>% mutate(tool=as.factor(tool))

tools <- c("iCluster", "lracluster", "mcca", "nemo", "pins", "snf",  "sumo")
names(tools) <- c('iClusterBayes', 'LRAcluster', 'MCCA', 'NEMO', 'PINSPlus', 'SNF', 'SUMO')
dataARI$tool <- unlist(lapply(dataARI$tool, function(x){tools[tools == x] %>% names()}))

dataARI$tool <- ordered(dataARI$tool, levels =c("NEMO", "SUMO", "LRAcluster", "MCCA", "PINSPlus", "SNF", "iClusterBayes"))

cairo_pdf("figureS2.pdf", width=8, height=5.5)
dataARI %>%
  ggplot() + 
  geom_ribbon(aes(x = var, ymin = minARI, ymax = maxARI, fill=tool), alpha=0.3) +
  geom_line(aes(x=var,y=medianARI, group=tool, color=tool), size=1) + 
  geom_point(aes(x=var,y=medianARI, group=tool, color=tool), size=1) + 
  facet_wrap(tool~., ncol=4) +
  theme_bw() + ylab("median ARI") + xlab("standard deviation") + 
  theme(legend.position = c(0.85, 0.25), legend.direction = "vertical", panel.grid.minor.y = element_blank(), 
        legend.title = element_blank(), text = element_text(size=12)) +
  scale_color_npg() + scale_fill_npg() + ylim(0,1)
dev.off()

dataARI %>%
  ggplot() + 
  geom_line(aes(x=var,y=medianARI, group=tool, color=tool), size=1.5) +
  geom_point(aes(x=var,y=medianARI, group=tool, color=tool), size=1.5) +
  theme_bw() + ylab("median ARI") + xlab("σ") + 
  theme(legend.position = "bottom", axis.text = element_text(size=12), strip.text.x = element_text(size = 12, face="bold"), 
        axis.title = element_text(size=12), legend.text=element_text(size=12), legend.title = element_blank()) +
  scale_color_npg()
