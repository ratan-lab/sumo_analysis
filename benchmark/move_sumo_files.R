
SUBTYPES.DATA = list(
  list(name='aml', only.primary=F, is.rna.seq=T, is.mirna.seq=T, display.name='AML'),
  list(name='breast', only.primary=T, is.rna.seq=T, is.mirna.seq=T, display.name='BIC'),
  list(name='colon', only.primary=T, is.rna.seq=T, is.mirna.seq=T, display.name='COAD'),
  list(name='gbm', only.primary=T, is.rna.seq=F, is.mirna.seq=F, display.name='GBM'),
  list(name='kidney', only.primary=T, is.rna.seq=T, is.mirna.seq=T, display.name='KIRC'),
  list(name='liver', only.primary=T, is.rna.seq=T, is.mirna.seq=T, display.name='LIHC'),
  list(name='lung', only.primary=T, is.rna.seq=T, is.mirna.seq=T, display.name='LUSC'),
  list(name='melanoma', only.primary=F, is.rna.seq=T, is.mirna.seq=T, display.name='SKCM'),
  list(name='ovarian', only.primary=T, is.rna.seq=T, is.mirna.seq=T, display.name='OV'),
  list(name='sarcoma', only.primary=T, is.rna.seq=T, is.mirna.seq=T, display.name='SARC'))

MAX.NUM.CLUSTERS = 15

results_dir <- "results"
sumo_files <- "sumo_files"

for (s in SUBTYPES.DATA){
  subtype <- s$name
  for (k in 2:MAX.NUM.CLUSTERS){
    #clustering
    fname <- file.path(sumo_files, subtype, paste0('k',k),'clusters.tsv')
    data <- read.table(fname, check.names = F, header = T, row.names = 1)
    clustering <- data$label
    names(clustering) <- rownames(data)
    outfile <- file.path(results_dir,paste(subtype, paste0('sumo',k), 'all', sep = "_"))
    save(clustering, file = outfile)
    #timing
    cmd <- paste('cp', file.path(results_dir,paste(subtype, 'sumo', 'all','timing', sep = "_")), file.path(results_dir,paste(subtype, paste0('sumo',k), 'all','timing', sep = "_")))
    cmd.return = system(cmd)
    stopifnot(cmd.return == 0)
  }
}

# ALGORITHM.NAMES = paste0('sumo', 2:MAX.NUM.CLUSTERS)
# ALGORITHM.DISPLAY.NAMES = as.list(paste0('SUMO', 2:MAX.NUM.CLUSTERS))
# names(ALGORITHM.DISPLAY.NAMES) = ALGORITHM.NAMES

