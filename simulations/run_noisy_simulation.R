source("noisy_simulation.R")

run.simulation <- function(){
  for (dir_path in c(SIMULATION.FILE.DIR, get.gauss.dir.path(), get.gamma.dir.path(), get.sumo.files.dir(),
                     get.clustering.results.dir.path("gamma"), get.clustering.results.dir.path("gauss"))){
    if (!dir.exists(dir_path)){
      dir.create(dir_path)
    }  
  }
  load.libraries()
  generate.original.dataset()
  
  GAMMA.FILES <- generate.gamma.data()
  SUMO.FILES.DIR <<- file.path(SIMULATION.FILE.DIR, "sumo_files", "gamma")
  results <- run.sampling(fnames=GAMMA.FILES, name="gamma")
  run.evaluation(results, outfile=file.path(SIMULATION.FILE.DIR,"gamma.tsv"))
  
  GAUSS.FILES <- generate.gauss.data()
  SUMO.FILES.DIR <<- file.path(SIMULATION.FILE.DIR, "sumo_files", "gauss")
  results <- run.sampling(fnames=GAUSS.FILES, name="gauss")
  run.evaluation(results, outfile=file.path(SIMULATION.FILE.DIR,"gauss.tsv"))
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
gamma_spl_k <- seq(0,5,0.25)
gamma_spl_theta <- 1

# gaussian_sampling
gauss_spl_mean <- 0
gauss_spl_std <- seq(0,3, 0.2)
gauss_spl_k <- 3
gauss_spl_theta <- 1

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

run.simulation()
