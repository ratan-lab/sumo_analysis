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
fraq = 0.2
base_mean = 0
base_sd =1
over_mean = c(2, 3) #first cluster have a fraction of features with mean=2 and sd=3 (set below) added to base distribution
over_sd = c(3, 2) 

# gamma_sampling
gamma_spl_mean <-2
gamma_spl_std <- 5
gamma_spl_k <- 1.5
gamma_spl_theta <- seq(0,10,0.5)

# gaussian_sampling
gauss_spl_mean <- seq(0,5, 0.25)
gauss_spl_std <- 5
gauss_spl_k <- 1.5
gauss_spl_theta <- 5

# simulation params
RANDOM.SEED <- 42
SIMULATION.FILE.DIR <- "noisy"

# benchmark params
MKL.BINARY.PATH ="run_MKL_DR/application"
MKL.ARGS.PATH = "MKL"
MCR.ROOT = "~/MATLAB/MATLAB_Runtime/v90"
MC.CORES <- 1
SUMO.PATH <- "sumo"
VARS.FNAME <- "../benchmark/set_vars.sh"
ALGORITHM.NAMES = c('snf','spectral', 'lracluster', 'pins', 'mcca', 'nemo', 'sumo', 'sumo_spectral') #'mkl'
ALGORITHM.DISPLAY.NAMES = as.list(c('SNF','Spectral', 'LRAcluster', 'PINS', 'MCCA', 'NEMO', 'SUMO', "SUMOspectral")) #'rMKL-LPP'
names(ALGORITHM.DISPLAY.NAMES) = ALGORITHM.NAMES
SUMO.FILES.DIR <- file.path(SIMULATION.FILE.DIR, "sumo_files")

run.simulation()
    
