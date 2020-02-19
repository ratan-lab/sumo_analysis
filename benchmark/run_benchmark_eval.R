source('nemo_benchmark/benchmark.R')
source('nemo_benchmark/NEMO.R')
source('SUMO.R')
source("benchmark.R")

# SUMO params
SUMO.FILES.DIR <- "sumo_files"
SUMO.PATH <- "sumo"
MAX.NUM.CLUSTERS <- 15

# benchmark params
DATASETS.PATH = "data"
RESULTS.DIR.PATH = "results"
CLINICAL.PARAMS.DIR = "data/clinical"
MC.CORES <- 25
VARS.FNAME <- "set_vars.sh"
RANDOM.SEED <- 42
OMIC.SUBSETS = list('multi_omics')
names(OMIC.SUBSETS) = c('all')

# run all tools
load.libraries()

repetitions <- 10
for (i in 1:repetitions){
  RESULTS.DIR.PATH = paste("results", i, sep = "_")
  SUMO.FILES.DIR <- paste("sumo_files", i, sep="_")
  RANDOM.SEED <- i
  
  ALGORITHM.NAMES = c('snf', 'lracluster', 'pins', 'mcca', 'nemo', 'sumo') #sumo
  ALGORITHM.DISPLAY.NAMES = as.list(c('SNF', 'LRAcluster', 'PINS', 'MCCA', 'NEMO', 'SUMO')) #SUMO
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
