source('nemo_benchmark/benchmark.R')
source('nemo_benchmark/NEMO.R')
source('SUMO.R')
source("benchmark.R")

# SUMO params
SUMO.PATH <- "sumo"

# benchmark params
DATASETS.PATH = "data"
RESULTS.DIR.PATH.MAIN = "results_eval"
CLINICAL.PARAMS.DIR = "data/clinical"
MC.CORES <- 20
VARS.FNAME <- "set_vars.sh"
RANDOM.SEED <- 42
OMIC.SUBSETS = list('multi_omics')
names(OMIC.SUBSETS) = c('all')

SUBTYPES.DATA = list(list(name='gbm', only.primary=T, is.rna.seq=F, is.mirna.seq=F, display.name='GBM'))
MAX.NUM.CLUSTERS <- 15

# run all tools
load.libraries()
dir.create(RESULTS.DIR.PATH.MAIN)

repetitions <- 10
for (i in 1:repetitions){
  RESULTS.DIR.PATH = file.path(RESULTS.DIR.PATH.MAIN, paste("results", i, sep = "_"))
  SUMO.FILES.DIR <- file.path(RESULTS.DIR.PATH.MAIN, paste("sumo_files", i, sep="_"))
  RANDOM.SEED <- i
  
  ALGORITHM.NAMES = c('snf', 'lracluster', 'pins', 'mcca', 'nemo', 'iCluster', 'sumo', 'cimlr')
  ALGORITHM.DISPLAY.NAMES = as.list(c('SNF', 'LRAcluster', 'PINS', 'MCCA', 'NEMO', 'iClusterBayes', 'SUMO', 'CIMLR'))
  names(ALGORITHM.DISPLAY.NAMES) = ALGORITHM.NAMES
  run.benchmark()
  
  # move sumo files
  if ('sumo' %in% ALGORITHM.NAMES){
    move.sumo.files()
    ALGORITHM.NAMES = c(ALGORITHM.NAMES, paste0('sumo', 2:MAX.NUM.CLUSTERS))
    ALGORITHM.DISPLAY.NAMES = as.list(c(unlist(ALGORITHM.DISPLAY.NAMES),paste0('SUMO', 2:MAX.NUM.CLUSTERS)))
    names(ALGORITHM.DISPLAY.NAMES) = ALGORITHM.NAMES
  }
  
  # calculate empirical survival and enriched clinical labels
  results <- analyze.benchmark()
  perform.all.analyses(results)
}
