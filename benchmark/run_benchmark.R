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
MC.CORES <- 20
VARS.FNAME <- "set_vars.sh"
RANDOM.SEED <- 42

ALGORITHM.NAMES = c('snf', 'lracluster', 'pins', 'mcca', 'nemo', 'iCluster', 'sumo', 'cimlr')
ALGORITHM.DISPLAY.NAMES = as.list(c('SNF', 'LRAcluster', 'PINS', 'MCCA', 'NEMO', 'iClusterBayes', 'SUMO', 'CIMLR'))
names(ALGORITHM.DISPLAY.NAMES) = ALGORITHM.NAMES

OMIC.SUBSETS = list('multi_omics')
names(OMIC.SUBSETS) = c('all')

# run other tools
load.libraries()
run.benchmark()
  
# move sumo files
if ('sumo' %in% ALGORITHM.NAMES){
  move.sumo.files()
  ALGORITHM.NAMES = c(ALGORITHM.NAMES,paste0('sumo', 2:MAX.NUM.CLUSTERS))
  ALGORITHM.DISPLAY.NAMES = as.list(c(unlist(ALGORITHM.DISPLAY.NAMES),paste0('SUMO', 2:MAX.NUM.CLUSTERS)))
  names(ALGORITHM.DISPLAY.NAMES) = ALGORITHM.NAMES
}
  
# calculate empirical survival and enriched clinical labels
results <- analyze.benchmark()
save(results, file.path(RESULTS.DIR.PATH, "benchmark_results_list.RData"))
perform.all.analyses(results)
print("DONE!")