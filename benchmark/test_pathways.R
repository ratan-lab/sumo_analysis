source('nemo_benchmark/benchmark.R')
source("benchmark.R")
library(tidyverse)
library(progeny)
library(broom)
library(ggpubr)

DATASETS.PATH = "data"
RESULTS.DIR.PATH = "results"

ALGORITHM.NAMES = c('snf', 'lracluster', 'pins', 'mcca', 'nemo', 'iCluster', 'sumo', 'cimlr')
ALGORITHM.DISPLAY.NAMES = as.list(c('SNF', 'LRAcluster', 'PINSPlus', 'MCCA', 'NEMO', 'iClusterBayes', 'SUMO', 'CIMLR'))
names(ALGORITHM.DISPLAY.NAMES) = ALGORITHM.NAMES

cancers <- c("AML", "BIC", "COAD", "GBM", "KIRC", "LIHC", "LUSC", "SKCM", "OV", "SARC")
selected <- c(12, 13, 7, 10, 6, 2, 9, 7, 9, 6)
names(selected) <- cancers

results <- NULL
for (subtype in SUBTYPES.DATA){
  # load expression data
  subtype.data <- get.raw.data(subtype.name = subtype$name,only.primary = subtype$only.primary, intersect.patients = T)
  gene_expr <- subtype.data[[1]] %>% as.matrix()
  print(subtype$display.name)
  print(dim(gene_expr))
  
  # fix gene names [annotate with HGNC symbols]
  fix.gene.names <- function(x){
    if(substring(x, 1, 3) == "X.." || substring(x, 1, 2) == "?|" ){
      return(x)
    } else {
      if (grepl("\\|", x)){
        return(strsplit(x,"\\|")[[1]][1])
      } 
      if (grepl("\\.", x)){
        return(strsplit(x,"\\.")[[1]][1])
      }
      return(x)
    }
  }
  gene_names <- sapply(rownames(gene_expr), fix.gene.names)
  rownames(gene_expr) <- gene_names
  
  pathways_all <- progeny(gene_expr, scale=TRUE, organism="Human", top=100, perm = 1)
  
  for (algorithm in ALGORITHM.NAMES){
    print(paste("testing", algorithm, "pathway enrichments using", subtype$name, "dataset"))

    # load sample labels
    fname <- file.path(RESULTS.DIR.PATH, paste(subtype$name, algorithm, "all", sep="_"))
    if (algorithm == "sumo"){
      fname <- file.path(RESULTS.DIR.PATH, paste(subtype$name, 
                                                 paste0(algorithm, selected[names(cancers) == subtype$display.name]), 
                                                 "all", sep="_"))
    }
    load(fname)
  
    labels <- clustering
    names(labels) <- colnames(gene_expr)
    pathways <- pathways_all[names(labels), ]
  
    # test enrichments
    data <- cbind(pathways, labels) %>% as_tibble(rownames = "sample") %>%
      gather(pathway, value, -sample, -labels) %>%
      mutate(labels = as.factor(labels))
  
    # data %>%
    #   ggboxplot(x="labels", y="value", facet.by="pathway") +
    #   stat_compare_means(label.y.npc = 0.02)
  
    test_assoc <- function(tbl) {
      res <- kruskal.test(tbl$value, tbl$labels)
      return(res$p.value)
    }
  
    enrichments <- data %>%
      group_by(pathway) %>%
      nest() %>%
      mutate(pval = map(data, test_assoc)) %>%
      select(-data) %>%
      unnest(cols = c(pval)) %>%
      mutate(tool=ALGORITHM.DISPLAY.NAMES[ALGORITHM.NAMES == algorithm][[1]],
             cancer=subtype$display.name)
  
    if (is.null(results)){
      results <- enrichments
    } else {
      results <- results %>% full_join(enrichments, by=c('pathway', 'pval', 'tool', 'cancer'))
    }
  }
}
total <- results %>% filter(pval < 0.05) %>% group_by(cancer, tool) %>% summarise(total=n())
write_tsv(results %>% spread(pathway, pval) %>% left_join(total) %>% arrange(cancer, tool), "benchmark_pathway_activity.tsv")
