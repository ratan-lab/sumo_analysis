source('nemo_benchmark/benchmark.R')
source('nemo_benchmark/NEMO.R')
source('SUMO.R')
source("benchmark.R")

# rMKL-LPP params
MKL.BINARY.PATH ="run_MKL_DR/application"
MKL.ARGS.PATH = "MKL"
MCR.ROOT = "~/MATLAB/MATLAB_Runtime/v90"

# # SUMO params
SUMO.FILES.DIR <- "sumo_files"
SUMO.PATH <- "sumo"

# benchmark params
DATASETS.PATH = "data"
RESULTS.DIR.PATH = "results"
CLINICAL.PARAMS.DIR = "data/clinical"
MC.CORES <- 25
VARS.FNAME <- "set_vars.sh"

ALGORITHM.NAMES = c('snf','spectral', 'lracluster', 'pins', 'mcca', 'nemo', 'mkl', 'sumo')
ALGORITHM.DISPLAY.NAMES = as.list(c('SNF','Spectral', 'LRAcluster', 'PINS', 'MCCA', 'NEMO', 'rMKL-LPP', 'SUMO'))
names(ALGORITHM.DISPLAY.NAMES) = ALGORITHM.NAMES

OMIC.SUBSETS = list('multi_omics')
names(OMIC.SUBSETS) = c('all')

# run all tools
load.libraries()
run.benchmark()

# calculate empirical survival and enriched clinical labels
results <- analyze.benchmark()
perform.all.analyses(results)

# plot results
p <- plot_results()
