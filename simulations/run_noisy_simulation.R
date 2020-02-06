source("noisy_simulation.R")

# base data info
nsamples =  200
nfeatures = 400
nclusters = 2
nlayers = 2
cluster_sd = 0.5

# simulation params
RANDOM.SEED <- 42
SIMULATION.FILE.DIR <- "noisy"

# benchmark params
# MKL.BINARY.PATH ="run_MKL_DR/application"
# MKL.ARGS.PATH = "MKL"
# MCR.ROOT = "~/MATLAB/MATLAB_Runtime/v90"
MC.CORES <- 4
SUMO.PATH <- "sumo"
GENERATE.DATA.SCRIPT <- "generate_data.py"
VARS.FNAME <- "../benchmark/set_vars.sh"
ALGORITHM.NAMES = c('snf','spectral', 'lracluster', 'pins', 'mcca', 'nemo', 'sumo', 'sumo_spectral') #'mkl'
ALGORITHM.DISPLAY.NAMES = as.list(c('SNF','Spectral', 'LRAcluster', 'PINS', 'MCCA', 'NEMO', 'SUMO', "SUMOspectral")) #'rMKL-LPP'
names(ALGORITHM.DISPLAY.NAMES) = ALGORITHM.NAMES
SUMO.FILES.DIR <- file.path(SIMULATION.FILE.DIR, "sumo_files")

###
prepare.simulation()

# constant layer
gauss_mean <-0
gauss_std <- 1.5
# sampled layer
gauss_spl_std <- seq(0,4, 2)
gauss_spl_mean <-0

SIM.FILES <- generate.two.gauss.layers(sampling="double_gauss")
results <- run.simulation(SIM.FILES, "double_gauss")
run.evaluation(results, outfile=file.path(SIMULATION.FILE.DIR, "double_gauss.tsv"))