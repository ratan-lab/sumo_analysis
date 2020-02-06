source("noisy_simulation.R")

prepare.simulation <- function(){
  for (dir_path in c(SIMULATION.FILE.DIR, get.sumo.files.dir())){
    if (!dir.exists(dir_path)){
      dir.create(dir_path)
    }  
  }
  load.libraries()
  generate.original.dataset()
}

run.simulation <- function(files, sampling_name){
  dir_name <- get.clustering.results.dir.path(sampling_name)
    if (!dir.exists(dir_name)){
      dir.create(dir_name)
    }  
  ORG.SUMO.FILES.DIR <- get.sumo.files.dir()
  SUMO.FILES.DIR <<- file.path(ORG.SUMO.FILES.DIR, sampling_name)
  results <- run.sampling(fnames=files, name=sampling_name)
  SUMO.FILES.DIR <<- ORG.SUMO.FILES.DIR
  return(results)
}


# original dataset info
nsamples =  200
nfeatures = 400
nclusters = 2
nlayers = 2
cluster_sd = 0.5

# gamma_sampling
gamma_spl_mean <-0
gamma_spl_std <- 1.5
gamma_spl_k <- seq(0,5,5)
gamma_spl_theta <- 1

# gaussian_sampling
gauss_spl_mean <- 0
gauss_spl_std <- seq(0,3,3)
gauss_spl_k <- 3
gauss_spl_theta <- 1

# # gaussian layers
# gamma_1_mean <-0
# gamma_1_std <- 1.5
# gamma_2_mean <-2
# gamma_2_std <- 1
# 
# # gaussian sampling
# gauss_spl_std <- seq(0,3, 0.2)

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
# ALGORITHM.NAMES = c('snf','spectral', 'lracluster', 'pins', 'mcca', 'nemo', 'sumo', 'sumo_spectral') #'mkl'
# ALGORITHM.DISPLAY.NAMES = as.list(c('SNF','Spectral', 'LRAcluster', 'PINS', 'MCCA', 'NEMO', 'SUMO', "SUMOspectral")) #'rMKL-LPP'
ALGORITHM.NAMES = c('snf','sumo', 'sumo_spectral')
ALGORITHM.DISPLAY.NAMES = as.list(c('SNF', 'SUMO', "SUMOspectral"))
names(ALGORITHM.DISPLAY.NAMES) = ALGORITHM.NAMES
SUMO.FILES.DIR <- file.path(SIMULATION.FILE.DIR, "sumo_files")

###
prepare.simulation()

GAMMA.FILES <- generate.gamma.data()
results.gamma <- run.simulation(GAMMA.FILES, "gamma")
run.evaluation(results.gamma, outfile=file.path(SIMULATION.FILE.DIR, "gamma.tsv"))

GAUSS.FILES <- generate.gauss.data()
results.gauss <- run.simulation(GAUSS.FILES, "gauss")
run.evaluation(results.gauss, outfile=file.path(SIMULATION.FILE.DIR, "gauss.tsv"))
