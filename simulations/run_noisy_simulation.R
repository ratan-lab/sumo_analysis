source("noisy_simulation.R")

# base data info
nsamples =  200
nfeatures = 400
nclusters = 2
nlayers = 2
cluster_sd = 0.5

# benchmark params
MC.CORES <- 4
SUMO.PATH <- "sumo"
GENERATE.DATA.SCRIPT <- "generate_data.py"
VARS.FNAME <- "../benchmark/set_vars.sh"
ALGORITHM.NAMES = c('snf', 'lracluster', 'pins', 'mcca', 'nemo', 'sumo')
ALGORITHM.DISPLAY.NAMES = as.list(c('SNF', 'LRAcluster', 'PINS', 'MCCA', 'NEMO', 'SUMO'))
names(ALGORITHM.DISPLAY.NAMES) = ALGORITHM.NAMES

# constant layer
gauss_mean <- 0
gauss_std <- 1.5
# sampled layer
gauss_spl_std <- seq(0,4, 0.2)
gauss_spl_mean <- 0

repetitions = 10
for (rep in 1:repetitions){
  RANDOM.SEED <- 42
  SIMULATION.FILE.DIR <- paste("noisy", rep, sep="_")
  SUMO.FILES.DIR <- file.path(SIMULATION.FILE.DIR, "sumo_files")
  prepare.simulation()
}

for (rep in 1:repetitions){
  print(paste0("#REP: ", rep))
  RANDOM.SEED <- rep
  SIMULATION.FILE.DIR <- paste("noisy", rep, sep="_")
  SUMO.FILES.DIR <- file.path(SIMULATION.FILE.DIR, "sumo_files")
  SIM.FILES <- generate.two.gauss.layers(sampling="double_gauss")
  results <- run.simulation(SIM.FILES, "double_gauss")
  run.evaluation(results, outfile=file.path(SIMULATION.FILE.DIR, "double_gauss.tsv"))
}
  