source('../benchmark/nemo_benchmark/benchmark.R')
source('../benchmark/benchmark.R')
source('../benchmark/nemo_benchmark/NEMO.R')
source('../benchmark/SUMO.R')

get.base.data.path <- function(){
  return(file.path(SIMULATION.FILE.DIR, "base_data.tsv"))
}

get.base.labels.path <- function(){
  return(file.path(SIMULATION.FILE.DIR, "base_labels.tsv"))
}

get.data.dir.path <- function(sampling){
  return(file.path(SIMULATION.FILE.DIR, paste(sampling, "data", sep="_")))
}

get.gamma.dir.path <- function(){
  return(file.path(SIMULATION.FILE.DIR, "gamma_data"))
}

get.gauss.dir.path <- function(){
  return(file.path(SIMULATION.FILE.DIR, "gauss_data"))
}

get.clustering.results.dir.path <- function(name){
  return(file.path(SIMULATION.FILE.DIR, paste0(name,"_results")))
}

get.generate.data.script <-function(){
  return(GENERATE.DATA.SCRIPT)
}

generate.original.dataset <- function(){
  # create base data and base labels
  cmd = paste("python3", get.generate.data.script(), "-rstate", RANDOM.SEED ,"-sd",cluster_sd, nsamples, nfeatures, nclusters, SIMULATION.FILE.DIR)
  cmd.return = system(cmd, intern = F)
  stopifnot(cmd.return == 0)
}

generate.gamma.data <- function(sampling="gamma"){
  base_data <- read.table(get.base.data.path())
  data_dir <- get.data.dir.path(sampling)
  if(!dir.exists(data_dir)){
    dir.create(data_dir)
  }
  set.seed(RANDOM.SEED)
  GAMMA.FILES <- list()
  
  #create midly noisy gauss layer
  gauss_fname <- file.path(data_dir, "gauss_layer.tsv")
  print(paste("Created:", gauss_fname))
  data <- base_data +  matrix(rnorm(n=nsamples*nfeatures, mean = gamma_spl_mean, sd = gamma_spl_std), ncol=nsamples, nrow=nfeatures)
  write.table(data, file=gauss_fname, sep = "\t", row.names = T, col.names = T)
  
  # create different gamma layers
  for (i in 1:length(gamma_spl_k)){
    fname <- file.path(get.gamma.dir.path(), paste0("gamma_k_", gamma_spl_k[i],".tsv"))
    print(paste("Created:", fname))
    data <- base_data + matrix(rgamma(n=nsamples*nfeatures, shape = gamma_spl_k[i], scale = gamma_spl_theta), ncol=nsamples, nrow=nfeatures)
    write.table(data, file=fname, sep = "\t", row.names = T, col.names = T)
    GAMMA.FILES <- append(GAMMA.FILES, list(list(layer1=gauss_fname, layer2=fname)))
  }
  
  return(GAMMA.FILES)
}

generate.gauss.data <- function(sampling="gauss"){
  base_data <- read.table(get.base.data.path())
  data_dir <- get.data.dir.path(sampling)
  if(!dir.exists(data_dir)){
    dir.create(data_dir)
  }
  set.seed(RANDOM.SEED)
  GAUSS.FILES <- list()
  
  #create midly noisy gamma layer
  gamma_fname <- file.path(get.data.dir.path(sampling), "gamma_layer.tsv")
  print(paste("Created:", gamma_fname))
  data <- base_data +  matrix(rgamma(n=nsamples*nfeatures, shape = gauss_spl_k, scale = gauss_spl_theta), ncol=nsamples, nrow=nfeatures)
  write.table(data, file=gamma_fname, sep = "\t", row.names = T, col.names = T)
  
  # create different gauss layers
  for (i in 1:length(gauss_spl_std)){
    fname <- file.path(get.gauss.dir.path(), paste0("gauss_std_", gauss_spl_std[i],".tsv"))
    print(paste("Created:", fname))
    data <- base_data + matrix(rnorm(n=nsamples*nfeatures, mean = gauss_spl_mean, sd = gauss_spl_std[i]), ncol=nsamples, nrow=nfeatures)
    write.table(data, file=fname, sep = "\t", row.names = T, col.names = T)
    GAUSS.FILES <- append(GAUSS.FILES, list(list(layer1=fname, layer2=gamma_fname)))
  }
  
  return(GAUSS.FILES)
}

get.raw.data <-function(data){
  layer1 <- read.table(data$layer1)
  attr(layer1, 'name') <- 'layer1'
  attr(layer1, 'fname') <- data$layer1
  layer2 <- read.table(data$layer2)
  attr(layer2, 'name') <- 'layer2'
  attr(layer2, 'fname') <- data$layer2
  return(list(layer1, layer2))
}

run.snf <- function(omics.list, subtype) {
  omics.list = lapply(omics.list, normalize.matrix)  
  alpha=0.5
  T.val=30
  num.neighbors = round(ncol(omics.list[[1]]) / 10)
  similarity.data = lapply(omics.list, function(x) {affinityMatrix(SNFtool::dist2(as.matrix(t(x)),as.matrix(t(x))), 
                                                                   num.neighbors, alpha)})
  W = SNF(similarity.data, num.neighbors, T.val)  
  num.clusters = nclusters
  clustering = spectralClustering(W, num.clusters)
  return(list(clustering=clustering))
}

run.spectral <- function(omics.list, subtype) {
  omics.list = lapply(omics.list, normalize.matrix)  
  concat.omics = do.call(rbind, omics.list)
  similarity.data = affinityMatrix(SNFtool::dist2(as.matrix(t(concat.omics)),as.matrix(t(concat.omics))), 20, 0.5)
  num.clusters = nclusters
  clustering = spectralClustering(similarity.data, num.clusters)
  return(list(clustering=clustering))
}

run.pins <- function(omics.list, subtype) {
  omics.transposed = lapply(omics.list, t)
  pins.ret = PINSPlus::SubtypingOmicsData(dataList=omics.transposed, kMax = nlayers)
  clustering = pins.ret$cluster1
  return(list(clustering=clustering))
}

run.mkl <- function(omics.list, subtype) {
  omics.list = lapply(omics.list, normalize.matrix)
  export.subtype.to.mkl(omics.list, subtype)
  bin.path = get.mkl.binary.path()
  subtype.dir = paste(get.mkl.arguments.path(), subtype, sep='/')
  paste(subtype.dir, 'kernels', sep='/')
  rundir <- getwd()
  setwd(bin.path)
  command = paste("./run_run_MKL_DR.sh", get.mcr.root.path(), paste(rundir, subtype.dir, 'kernels', sep='/'),
                  paste(rundir, subtype.dir, 'ids', sep='/'),
                  paste(rundir, subtype.dir, 'output', sep='/'),
                  '9', '5')
  command.return = system(command)
  setwd(rundir)
  stopifnot(command.return == 0)
  clustering = get.mkl.clustering(subtype)
  return(list(clustering=clustering))
}

run.mcca <- function(omics.list, subtype, known.num.clusters=nclusters, penalty=NULL, rep.omic=1) {# 2D noisy simulation
  omics.list = lapply(omics.list, normalize.matrix)  
  max.dim = min(nlayers, nrow(omics.list[[1]]))
  omics.transposed = lapply(omics.list, t)
  cca.ret = PMA::MultiCCA(omics.transposed, type="ordered", ncomponents = 1, penalty=penalty) #type="standard", ncomponents = max.dim
  sample.rep = omics.transposed[[rep.omic]] %*% cca.ret$ws[[rep.omic]]
  cca.clustering = kmeans(sample.rep, known.num.clusters, iter.max=100, nstart=30)$cluster  
  return(list(clustering=cca.clustering))
}

run.nemo <- function(omics.list, subtype.data, num.clusters=nclusters, is.missing.data=F, num.neighbors=NA) {# 2D noisy simulation
  omics.list = lapply(omics.list, normalize.matrix)  
  clustering = nemo.clustering(omics.list, num.clusters, num.neighbors)
  return(list(clustering=clustering))
}

run.lracluster <- function(omics.list, subtype) {
  dim.range = 1:nlayers
  all.clustering.results = list()
  
  omics.matrix.list = lapply(omics.list, as.matrix)
  dimension = nlayers# 2D noisy simulation
  print(paste('running lra cluster for dimension', dimension))
  data.names = c('gauss','gamma')
  clustering.results = LRAcluster(omics.matrix.list, 
                                  rep('gaussian', length(omics.list)), 
                                  dimension=dimension, data.names)
  solution = clustering.results$coordinate
  
  num.clusters = nclusters# 2D noisy simulation
  print(paste('running kmeans in lra cluster for num clusters', num.clusters))
  clustering = kmeans(t(solution), num.clusters, iter.max=100, nstart=60)$cluster
  return(list(clustering=clustering))
}

run.sumo <- function(omics.list, subtype, num.clusters=nclusters, mc.cores=get.mc.cores(), file_dir= get.sumo.files.dir()){
  if (!dir.exists(file_dir)){
    dir.create(file_dir)
  }
  prepare_files <- c()
  sampling <- strsplit(subtype, "_")[[1]][1]
  
  for (i in 1:length(omics.list)){
    if (attr(omics.list[[i]], "name") == sampling){
      fname <- file.path(get.sumo.files.dir(), paste0(subtype, "_std.tsv"))
      omic <- t(scale(t(as.matrix(omics.list[[i]]))))
      omics.list[[i]] <- omic
      write.table(omic, fname, sep="\t", row.names = T, col.names = T)
      prepare_files <- c(prepare_files, fname)
    } else {
      fname <- paste0(strsplit(attr(omics.list[[i]], "fname"), "\\.")[[1]][1], "_std.tsv" )
      if (!file.exists(fname)){
        omic <- t(scale(t(as.matrix(omics.list[[i]]))))
        omics.list[[i]] <- omic
        write.table(omic, fname, sep="\t", row.names = T, col.names = T)
      }
      prepare_files <- c(prepare_files, fname)
    }
  }
  
  outfile = file.path(file_dir, paste0("prepared.", subtype, ".npz"))
  sumo.prepare(omics.list, subtype.data, prepare_files, outfile=outfile)
  stopifnot(file.exists(outfile))
  
  outdir = file.path(file_dir, subtype)
  clustering = sumo.clustering(num.clusters, outfile, outdir, mc.cores)
  return(list(clustering=clustering))
}

run.sumo_spectral <- function(omics.list, subtype, num.clusters=nclusters, mc.cores=get.mc.cores(), file_dir= get.sumo.files.dir()){
  if (!dir.exists(file_dir)){
    dir.create(file_dir)
  }
  prepare_files <- c()
  sampling <- strsplit(subtype, "_")[[1]][1]
  
  for (i in 1:length(omics.list)){
    if (attr(omics.list[[i]], "name") == sampling){
      fname <- file.path(get.sumo.files.dir(), paste0(subtype, "_std.tsv"))
      omic <- t(scale(t(as.matrix(omics.list[[i]]))))
      omics.list[[i]] <- omic
      write.table(omic, fname, sep="\t", row.names = T, col.names = T)
      prepare_files <- c(prepare_files, fname)
    } else {
      fname <- paste0(strsplit(attr(omics.list[[i]], "fname"), "\\.")[[1]][1], "_std.tsv" )
      if (!file.exists(fname)){
        omic <- t(scale(t(as.matrix(omics.list[[i]]))))
        omics.list[[i]] <- omic
        write.table(omic, fname, sep="\t", row.names = T, col.names = T)
      }
      prepare_files <- c(prepare_files, fname)
    }
  }
  
  outfile = file.path(file_dir, paste0("prepared.", subtype, ".npz"))
  sumo.prepare(omics.list, subtype.data, prepare_files, outfile=outfile)
  stopifnot(file.exists(outfile))
  
  outdir = file.path(file_dir, subtype)
  clustering = sumo.clustering_spectral(num.clusters, outfile, outdir, mc.cores)
  return(list(clustering=clustering))
}

sumo.clustering <- function(k, fname, outdir, mc.cores = get.mc.cores()){
  log <- paste0(paste(strsplit(fname, "\\.")[[1]][1:2], collapse = '.'), ".log")
  command <- paste(get.sumo.path(), "run","-t", mc.cores,"-log DEBUG", "-logfile", log, fname, k, outdir)
  print(command)
  command.return = system(command)
  stopifnot(command.return == 0)
  
  #TODO select the best results and point to them below
  results_file <- file.path(outdir, paste0("k",k), "sumo_results.npz")
  np <- import("numpy")
  npz <- np$load(results_file, allow_pickle = T)
  clusters <- npz$f[["clusters"]]
  
  clustering <- unlist(clusters[,2])
  names(clustering) <- unlist(clusters[,1])
  return(clustering)
}


sumo.clustering_spectral <- function(k, fname, outdir, mc.cores = get.mc.cores()){
  log <- paste0(paste(strsplit(fname, "\\.")[[1]][1:2], collapse = '.'), ".log")
  command <- paste(get.sumo.path(), "run","-method","spectral","-t", mc.cores,"-log DEBUG", "-logfile", log, fname, k, outdir)
  print(command)
  command.return = system(command)
  stopifnot(command.return == 0)
  
  #TODO select the best results and point to them below
  results_file <- file.path(outdir, paste0("k",k), "sumo_results.npz")
  np <- import("numpy")
  npz <- np$load(results_file, allow_pickle = T)
  clusters <- npz$f[["clusters"]]
  
  clustering <- unlist(clusters[,2])
  names(clustering) <- unlist(clusters[,1])
  return(clustering)
}


run.sampling <- function(fnames, name) {
  result_files <- list()
  for (i in 1:length(fnames)) {
    current.data = fnames[[i]]
    subtype = paste(name, i, sep="_")
    subtype.raw.data = get.raw.data(current.data)
    stopifnot(all(colnames(subtype.raw.data[[1]]) == colnames(subtype.raw.data[[2]])))
    samples <- colnames(subtype.raw.data[[1]])
    for (algorithm.name in ALGORITHM.NAMES) {
        set.seed(RANDOM.SEED)
        print(paste('data', subtype, 'running algorithm', algorithm.name))
        clustering.path = file.path(get.clustering.results.dir.path(name), paste0(subtype, "_", algorithm.name, ".tsv"))
        if (!file.exists(clustering.path)) {
          algorithm.func.name = paste0('run.', algorithm.name)
          algorithm.func = get(algorithm.func.name)
          algorithm.ret = algorithm.func(subtype.raw.data, subtype)
          clustering = algorithm.ret$clustering
          print('before saving')
          if(!is.null(names(clustering))){
            write.table(data.frame(sample=names(clustering), label=clustering), file=clustering.path, sep="\t", row.names = F, col.names = T)
          } else {
            # assume the same sample order
            write.table(data.frame(sample=samples, label=clustering), file=clustering.path, sep="\t", row.names = F, col.names = T)
          }
          result_files <- append(result_files, list(list(fname=clustering.path, subtype=subtype, algorithm=algorithm.name)))
        }
    }
  }
  return(result_files)
}


run.evaluation <- function(results, outfile){
  algorithms <- c()
  subtypes <- c()
  vals <- c()
  metrics <- c()
  for (r in results){
     stopifnot(file.exists(r$fname))
    cmd <- paste(get.sumo.path(), 'evaluate', r$fname, get.base.labels.path())
    cmd.return <- system(cmd, intern=T)
    for (metric in c('ARI', 'NMI','purity')){
      line <- cmd.return[grepl(metric, cmd.return)]
      val <- strsplit(line,'\t')[[1]][2]
      algorithms <- c(algorithms, r$algorithm)
      subtypes <- c(subtypes, r$subtype)
      vals <- c(vals, val)
      metrics <- c(metrics, metric)
    }
  }
  data <- data.frame(algoritm=algorithms, subtype=subtypes, metric = metrics, val = vals)
  write.table(data, file=outfile, sep = "\t", row.names = F, col.names = T)
}
