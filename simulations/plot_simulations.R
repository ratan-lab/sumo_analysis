# Create Fig1 (noisy simulation and missing simulation) & FigS2 (noisy simulation detailed)
library(tidyverse)
library(ggsci)
library(cowplot)

noisy_results_fname <- "noisy/double_gauss.tsv" 
missing_results_fname <- "missing/missing.tsv"
gauss_spl_std <- seq(0,4, 0.2)

# create color palettes
ALGORITHM.NAMES <- c("nemo", "sumo", "lracluster", "mcca", "pins", "snf", "iCluster", "cimlr")
names(ALGORITHM.NAMES) <- c("NEMO", "SUMO","LRAcluster", "MCCA", "PINSPlus", "SNF", "iClusterBayes", "CIMLR")
tool_pal_color <- pal_npg("nrc")(10)[c(1:7,10)]
names(tool_pal_color) <- names(ALGORITHM.NAMES)
tool_pal_fill <- pal_npg("nrc", alpha = 0.7)(10)[c(1:7,10)]
names(tool_pal_fill) <- names(ALGORITHM.NAMES)

# load noisy simulation data
map_vals <- function(val){
  RANGE <- gauss_spl_std
  return(RANGE[as.integer(val)])
}

noisy_data <- read_tsv(noisy_results_fname, col_types = cols()) 

noisy_data <- noisy_data %>% 
  separate(subtype, c("sampling", "var"), sep = "_(?=[:digit:])") %>%
  mutate(var = map_vals(var)) %>%
  filter(metric == "ARI")

noisy_data <- noisy_data %>% 
  group_by(tool, var) %>% 
  summarise(medianARI= median(val), minARI=min(val), maxARI=max(val))  %>% 
  right_join(noisy_data, by = c("tool", "var")) %>%
  ungroup() %>%
  mutate(tool=as.factor(tool))

noisy_data$tool <- unlist(lapply(noisy_data$tool, function(x){ALGORITHM.NAMES[ALGORITHM.NAMES == x] %>% names()}))
noisy_data$tool <- ordered(noisy_data$tool, levels = names(ALGORITHM.NAMES))

# load missing data
missing_data <- read_tsv(missing_results_fname, col_types = cols())

missing_data <- missing_data %>% 
  filter(metric == "ARI") %>%
  mutate(tool = as.factor(tool))

levels(missing_data$tool) <- unlist(lapply(levels(missing_data$tool), function(x){ALGORITHM.NAMES[ALGORITHM.NAMES == x] %>% names()}))

missing_data <- missing_data%>%
  group_by(tool, metric, fraction) %>%
  summarize(median=median(val), min=min(val), max=max(val)) %>%
  right_join(missing_data, by = c("tool", "metric", "fraction")) %>%
  mutate(fraction=fraction*100)

# create figures

cairo_pdf("figureS2.pdf", width=8, height=5.5)

  noisy_data %>%
    ggplot() + 
    geom_ribbon(aes(x = var, ymin = minARI, ymax = maxARI, fill=tool), alpha=0.3) +
    geom_line(aes(x=var,y=medianARI, group=tool, color=tool), size=1) + 
    geom_point(aes(x=var,y=medianARI, group=tool, color=tool), size=1) + 
    facet_wrap(tool~., ncol=4) +
    theme_bw() + ylab("median ARI") + xlab("standard deviation") + 
    theme(legend.position=c(0.85,0.25), 
          legend.direction = "vertical", 
          panel.grid.minor.y = element_blank(), 
          legend.title = element_blank(), text = element_text(size=12)) +
    scale_color_manual(values=tool_pal_color) + scale_fill_manual(values=tool_pal_fill) + ylim(0,1)

dev.off()

cairo_pdf("figure1.pdf", width = 11.5, height = 5.5) 

  p1 <- noisy_data %>%
    ggplot() + 
    geom_line(aes(x=var,y=medianARI, group=tool, color=tool), size=1.5) +
    geom_point(aes(x=var,y=medianARI, group=tool, color=tool), size=1.5) +
    theme_bw() + ylab("median ARI") + xlab("standard deviation") + 
    theme(legend.position = "bottom", axis.text = element_text(size=12), strip.text.x = element_text(size = 12, face="bold"), 
          axis.title = element_text(size=12), legend.text=element_text(size=12), legend.title = element_blank()) +
    scale_color_manual(values=tool_pal_color) 
  
  p2 <- missing_data %>%
    ggplot() + 
    geom_line(aes(x=as.factor(fraction), y=median, group=tool, color=tool), size=1.5) +
    geom_point(aes(x=as.factor(fraction),y=median, group=tool, color=tool), size=1.5) + 
    geom_ribbon(aes(x=as.factor(fraction), ymin=min, ymax=max, group=tool, fill=tool), alpha=0.3, show.legend = F) +
    theme_bw() + theme(legend.position = "bottom", axis.text = element_text(size=14), axis.title = element_text(size=14),
                       legend.text=element_text(size=14), strip.text.x = element_text(size = 14, face="bold")) +
    ylab("median ARI") + xlab("percent of missing samples") + labs(color="") + scale_color_npg() + scale_fill_npg()
  
  legend <- get_legend(p1 + guides(color = guide_legend(nrow = 1)) + 
                         theme(legend.position = "bottom"))
  
  prow <- plot_grid(p1 + theme(legend.position="none"),
                    p2 + theme(legend.position="none"),
                    labels = c('A','B'))
  
  plot_grid(prow, legend, ncol = 1, rel_heights = c(1, .1))

dev.off()
