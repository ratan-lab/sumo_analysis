source("missing_simulation.R")

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
ALGORITHM.NAMES = c( 'nemo', 'sumo')
ALGORITHM.DISPLAY.NAMES = as.list(c('NEMO', 'SUMO'))
names(ALGORITHM.DISPLAY.NAMES) = ALGORITHM.NAMES

# noise in layers
gauss_mean_1 <- 1
gauss_std_1 <- 1
gauss_mean_2 <- 0
gauss_std_2 <- 1.5

repetitions = 10
fraq <- seq(0,0.9,0.1) #fraction of samples removed from first layer

all_results <- list()
for (rep in 1:repetitions){
  # prepare simulation
  RANDOM.SEED <- 42
  SIMULATION.FILE.DIR <- paste("missing", rep, sep="_")
  prepare.simulation()
  
  # noisy data
  RANDOM.SEED <- rep
  layer1 <- file.path(SIMULATION.FILE.DIR, "layer1.tsv")
  layer2 <- file.path(SIMULATION.FILE.DIR, "layer2.tsv")
  generate.noisy.data(gauss_mean_1, gauss_std_1, outfile=layer1)
  generate.noisy.data(gauss_mean_2, gauss_std_2, outfile=layer2)
  
  for (fr in fraq){
    if (fr >= 1 | fr < 0){
      stop(paste("Cannot remove",fr*100,"% of samples!"))
    }
    # subsample data
    print(paste0("#REP: ", rep, " (", fr*100, "% samples removed from first layer)"))
    SUMO.FILES.DIR <- file.path(SIMULATION.FILE.DIR, paste0("sumo_files_", fr))
    SIM.FILES <- generate.missing.data(fraction=fr, layer1=layer1, layer2=layer2)
    # run simulation
    results <- run.simulation(SIM.FILES, sampling=paste0(fr))
    all_results <- append(all_results, results)
  }
}

run.evaluation(all_results, outfile="missing.tsv")
print('DONE')