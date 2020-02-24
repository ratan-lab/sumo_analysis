source("subsampling_simulation.R")

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

# noise in layers
gauss_mean_1 <- 0
gauss_std_1 <- 1.5
gauss_mean_2 <- 1
gauss_std_2 <- 1

repetitions = 100
fraq <- seq(0,0.9,0.1)

# prepare simulation 
SIMULATION.FILE.DIR <- "subsampling"
RANDOM.SEED <- 42 # used only for base data creation (not for sampling)
prepare.simulation(rseed = RANDOM.SEED)
LAYER1 <- file.path(SIMULATION.FILE.DIR, "layer1.tsv")
LAYER2 <- file.path(SIMULATION.FILE.DIR, "layer2.tsv")
generate.noisy.data(gauss_mean_1, gauss_std_1, outfile=LAYER1)
generate.noisy.data(gauss_mean_2, gauss_std_2, outfile=LAYER2)

# run simulation
all_results <- list()
for (rep in 1:repetitions){
  SIMULATION.FILE.DIR <- file.path("subsampling", paste("subsampling", rep, sep="_"))
  for (fr in fraq){
    if (fr >= 1 | fr < 0){
      stop(paste("Cannot remove",fr*100,"% of samples!"))
    }
    # subsample data
    print(paste0("#REP: ", rep, " (", fr*100, "% samples removed from both layers)"))
    SUMO.FILES.DIR <- file.path(SIMULATION.FILE.DIR, paste0("sumo_files_", fr))
    SIM.FILES <- generate.subsampled.data(fraction=fr, layer1=LAYER1, layer2=LAYER2, outdir=SIMULATION.FILE.DIR)
    # run tools
    results <- run.simulation(SIM.FILES, sampling=paste0(fr))
    all_results <- append(all_results, results)
  }
}

SIMULATION.FILE.DIR <- "subsampling"
run.evaluation(all_results, outfile=file.path("subsampling","subsampling.tsv"))
