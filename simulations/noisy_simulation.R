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

generate.two.gauss.layers <- function(sampling){
  base_data <- read.table(get.base.data.path())
  data_dir <- get.data.dir.path(sampling)
  if(!dir.exists(data_dir)){
    dir.create(data_dir)
  }
  set.seed(RANDOM.SEED)
  FILE.LIST <- list()
  
  #create constant layer
  constant_fname <- file.path(data_dir, "constant_layer.tsv")
  print(paste("Created:", constant_fname))
  data <- base_data +  matrix(rnorm(n=nsamples*nfeatures, mean = gauss_mean, sd = gauss_std), ncol=nsamples, nrow=nfeatures)
  write.table(data, file=constant_fname, sep = "\t", row.names = T, col.names = T)
  
  # create different gauss layers
  for (i in 1:length(gauss_spl_std)){
    fname <- file.path(data_dir, paste0("gauss_std_", gauss_spl_std[i],".tsv"))
    print(paste("Created:", fname))
    data <- base_data + matrix(rnorm(n=nsamples*nfeatures, mean = gauss_spl_mean, sd = gauss_spl_std[i]), ncol=nsamples, nrow=nfeatures)
    write.table(data, file=fname, sep = "\t", row.names = T, col.names = T)
    FILE.LIST <- append(FILE.LIST, list(list(layer1=fname, layer2=constant_fname)))
  }
  
  return(FILE.LIST)
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
  start = Sys.time()
  omics.list = lapply(omics.list, normalize.matrix)  
  alpha=0.5
  T.val=30
  num.neighbors = round(ncol(omics.list[[1]]) / 10)
  similarity.data = lapply(omics.list, function(x) {affinityMatrix(SNFtool::dist2(as.matrix(t(x)),as.matrix(t(x))), 
                                                                   num.neighbors, alpha)})
  W = SNF(similarity.data, num.neighbors, T.val)  
  num.clusters = nclusters
  clustering = spectralClustering(W, num.clusters)
  time.taken = as.numeric(Sys.time() - start, units='secs')
  return(list(clustering=clustering, timing=time.taken))
}

run.pins <- function(omics.list, subtype) {
  start = Sys.time()
  omics.transposed = lapply(omics.list, t)
  pins.ret = PINSPlus::SubtypingOmicsData(dataList=omics.transposed, kMax = nlayers)
  clustering = pins.ret$cluster1
  time.taken = as.numeric(Sys.time() - start, units='secs')
  return(list(clustering=clustering, timing=time.taken))
}

run.mcca <- function(omics.list, subtype, known.num.clusters=nclusters, penalty=NULL, rep.omic=1) {# 2D noisy simulation
  start = Sys.time()
  omics.list = lapply(omics.list, normalize.matrix)  
  max.dim = min(nlayers, nrow(omics.list[[1]]))
  omics.transposed = lapply(omics.list, t)
  cca.ret = PMA::MultiCCA(omics.transposed, type="ordered", ncomponents = 1, penalty=penalty) #type="standard", ncomponents = max.dim
  sample.rep = omics.transposed[[rep.omic]] %*% cca.ret$ws[[rep.omic]]
  cca.clustering = kmeans(sample.rep, known.num.clusters, iter.max=100, nstart=30)$cluster  
  time.taken = as.numeric(Sys.time() - start, units='secs')
  return(list(clustering=cca.clustering, timing=time.taken))
}

run.nemo <- function(omics.list, subtype.data, num.clusters=nclusters, is.missing.data=F, num.neighbors=NA) {# 2D noisy simulation
  start = Sys.time()
  omics.list = lapply(omics.list, normalize.matrix)  
  clustering = nemo.clustering(omics.list, num.clusters, num.neighbors)
  time.taken = as.numeric(Sys.time() - start, units='secs')
  return(list(clustering=clustering, timing=time.taken))
}

run.lracluster <- function(omics.list, subtype) {
  start = Sys.time()
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
  time.taken = as.numeric(Sys.time() - start, units='secs')
  return(list(clustering=clustering, timing=time.taken))
}

run.sumo <- function(omics.list, subtype, num.clusters=nclusters, mc.cores=get.mc.cores(), file_dir= get.sumo.files.dir()){
  start = Sys.time()
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
  time.taken = as.numeric(Sys.time() - start, units='secs')
  return(list(clustering=clustering, timing=time.taken))
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
        timing.path = file.path(get.clustering.results.dir.path(name), paste0(subtype, "_",algorithm.name, '_timing'))
        if (!file.exists(clustering.path)) {
          algorithm.func.name = paste0('run.', algorithm.name)
          algorithm.func = get(algorithm.func.name)
          algorithm.ret = algorithm.func(subtype.raw.data, subtype)
          clustering = algorithm.ret$clustering
          timing = algorithm.ret$timing
          print('before saving')
          save(timing, file = timing.path)
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
  data <- data.frame(tool=algorithms, subtype=subtypes, metric = metrics, val = vals)
  write.table(data, file=outfile, sep = "\t", row.names = F, col.names = T)
}

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
