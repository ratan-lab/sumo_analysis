library(tidyverse)
library(ggsci)

data <- read_tsv('subsampling.tsv')
dataARI <- data %>% filter(metric == "ARI") 

dataARI%>%
  group_by(tool, metric, subtype) %>%
  summarize(median=median(val), min=min(val), max=max(val)) %>%
  full_join(dataARI) %>%
  mutate(subtype=subtype*100) %>%
  ggplot() + geom_line(aes(x=as.factor(subtype), y=median, group=tool, color=tool), size=1.5) + 
  # geom_errorbar(aes(x=as.factor(subtype),y=median, ymin=min, ymax=max, group=tool,color=tool), alpha=0.05, size= 1) +
  theme_bw() + theme(legend.position = "bottom", axis.text = element_text(size=12), axis.title = element_text(size=12),
                     legend.text=element_text(size=12), strip.text.x = element_text(size = 12, face="bold")) + 
  ylab("median ARI") + xlab("percent of samples removed from both layers") + labs(color="") + scale_color_npg()

dataARI%>%
  group_by(tool, metric, subtype) %>%
  summarize(median=median(val), min=min(val), max=max(val)) %>%
  full_join(dataARI) %>%
  mutate(subtype=subtype*100) %>%
  ggplot() + geom_line(aes(x=as.factor(subtype), y=median, group=tool, color=tool), size=1) +
  geom_ribbon(aes(x=as.factor(subtype), ymin=min, ymax=max, group=tool, fill=tool), alpha=0.3, show.legend = F) +
  facet_wrap(tool~.) +
  theme_bw() + theme(legend.position = "bottom", axis.text = element_text(size=12), axis.title = element_text(size=12),
                     legend.text=element_text(size=12), strip.text.x = element_text(size = 12, face="bold")) +
  ylab("median ARI") + xlab("percent of samples removed from both layers") + labs(color="") + scale_color_npg() + scale_fill_npg()

####
library(Rtsne)
library(ggpubr)

data1 <- read.table("subsampling_1/layer1.tsv")
data2 <- read.table("subsampling_1/layer2.tsv")
labels <- read_tsv("subsampling_1/base_labels.tsv")
stopifnot(all(colnames(data1) == labels$sample) & all(colnames(data1) == colnames(data2)))

# first layer 
res <- Rtsne(t(as.matrix(data1)))
tsne_plot <- data.frame(x = res$Y[,1], y = res$Y[,2], col = labels$label)
ggplot(tsne_plot) + geom_point(aes(x=x, y=y, fill=as.factor(col)), size=4, pch=21,colour="black") + 
  xlab("tSNE1") + ylab("tSNE2") + labs(color="clusters") + theme_bw() + 
  theme(legend.position = "null", axis.text = element_text(size=12), axis.title = element_text(size=12),
        legend.text=element_text(size=12), strip.text.x = element_text(size = 12, face="bold"))
p1 <- last_plot()

# second layer
res <- Rtsne(t(as.matrix(data2)))
tsne_plot <- data.frame(x = res$Y[,1], y = res$Y[,2], col = labels$label)
ggplot(tsne_plot) + geom_point(aes(x=x, y=y, fill=as.factor(col)), size=4, pch=21,colour="black") + 
  xlab("tSNE1") + ylab("tSNE2") + labs(color="clusters") + theme_bw() + 
  theme(legend.position = "null", axis.text = element_text(size=12), axis.title = element_text(size=12),
        legend.text=element_text(size=12), strip.text.x = element_text(size = 12, face="bold"))
p2 <- last_plot()

ggarrange(p1,p2)
