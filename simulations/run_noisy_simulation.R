source("noisy_simulation.R")

# base data info
nsamples =  200
nfeatures = 400
nclusters = 2
nlayers = 2
cluster_sd = 0.5

# benchmark params
MC.CORES <- 15
SUMO.PATH <- "sumo"
GENERATE.DATA.SCRIPT <- "generate_data.py"
VARS.FNAME <- "../benchmark/set_vars.sh"
ALGORITHM.NAMES = c('snf', 'lracluster', 'pins', 'mcca', 'nemo', 'iCluster', 'sumo', 'cimlr')
ALGORITHM.DISPLAY.NAMES = as.list(c('SNF', 'LRAcluster', 'PINS', 'MCCA', 'NEMO', 'iClusterBayes', 'SUMO', 'CIMLR'))
names(ALGORITHM.DISPLAY.NAMES) = ALGORITHM.NAMES

# constant layer
gauss_mean <- 0
gauss_std <- 1.5
# sampled layer
gauss_spl_std <- seq(0,4, 0.2)
gauss_spl_mean <- 0

RANDOM.SEED <- 42 # used only for base data creation (not for noise generation)
repetitions = 100

# prepare simulation
MAIN.SIMULATION.FILE.DIR <- "noisy"
prepare.simulation(RANDOM.SEED)

# run simulation
all_results <- list()
for (rep in 1:repetitions){
  print(paste0("#REP: ", rep))
  SIMULATION.FILE.DIR <- file.path(MAIN.SIMULATION.FILE.DIR, paste(MAIN.SIMULATION.FILE.DIR, rep, sep="_"))
  if (dir.exists(SIMULATION.FILE.DIR)){
    print(paste("Directory", SIMULATION.FILE.DIR, "already exists!"))
  } else {
    dir.create(SIMULATION.FILE.DIR)
  }
  SUMO.FILES.DIR <- file.path(SIMULATION.FILE.DIR, "sumo_files")
  # SIM.FILES <- generate.two.gauss.layers(sampling="double_gauss")
  SIM.FILES <- get.two.gauss.layers(sampling="double_gauss") # use when re-running the simulation on existing files
  results <- run.simulation(SIM.FILES, "double_gauss")
  all_results <- append(all_results, results)
}

run.evaluation(all_results, outfile=file.path(MAIN.SIMULATION.FILE.DIR, "double_gauss.tsv"))
