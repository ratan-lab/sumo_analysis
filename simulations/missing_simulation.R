source('../benchmark/nemo_benchmark/benchmark.R')
source('../benchmark/benchmark.R')
source('../benchmark/nemo_benchmark/NEMO.R')
source('../benchmark/SUMO.R')

get.generate.data.script <-function(){
  return(GENERATE.DATA.SCRIPT)
}

get.base.data.path <- function(){
  return(file.path(SIMULATION.FILE.DIR, "base_data.tsv"))
}

get.base.labels.path <- function(){
  return(file.path(SIMULATION.FILE.DIR, "base_labels.tsv"))
}

get.clustering.results.dir.path <- function(){
  return(file.path(SIMULATION.FILE.DIR, paste0("results")))
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


run.iCluster <- function(omics.list, subtype.data) {
  start = Sys.time()
  dev.ratios = c()
  icluster.rets = list()
  icluster.ret = iClusterPlus::tune.iClusterBayes(cpus=get.mc.cores(), t(omics.list[[1]]), t(omics.list[[2]]), 
                                                  K=nclusters, type=rep('gaussian', 2))$fit
  dev.ratios = icluster.ret[[1]]$dev.ratio
  
  print('dev.ratios are:')
  print(dev.ratios)
  
  optimal.solution = icluster.ret[[which.max(dev.ratios)]]
  time.taken = as.numeric(Sys.time() - start, units='secs')
  return(list(clustering=optimal.solution$clusters, 
              timing=time.taken))
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
  
  for (i in 1:length(omics.list)){
    data <- omics.list[[i]]
    name <- attr(data, "name")
    fname <- file.path(get.sumo.files.dir(), paste0(name, "_std.tsv"))
    
    omic <- t(scale(t(as.matrix(data))))
    omics.list[[i]] <- omic
    write.table(omic, fname, sep="\t", row.names = T, col.names = T)
    prepare_files <- c(prepare_files, fname)
  }
  
  outfile = file.path(file_dir, paste0("prepared.", subtype, ".npz"))
  sumo.prepare(omics.list, subtype, prepare_files, outfile=outfile)
  stopifnot(file.exists(outfile))
  
  outdir = file.path(file_dir, paste0("run_", subtype))
  clustering = sumo.clustering(num.clusters, outfile, outdir, mc.cores)
  time.taken = as.numeric(Sys.time() - start, units='secs')
  return(list(clustering=clustering, timing=time.taken))
}

sumo.clustering <- function(k, fname, outdir, mc.cores = get.mc.cores()){
  log <- paste0(outdir, ".log")
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

prepare.simulation <- function(rseed){
  if (!dir.exists(SIMULATION.FILE.DIR)){
    dir.create(SIMULATION.FILE.DIR)
  }
  load.libraries()
  generate.original.dataset(rseed)
}

generate.original.dataset <- function(rseed){
  # create base data and base labels
  cmd = paste("python3", get.generate.data.script(), "-rstate", rseed ,"-sd",cluster_sd, nsamples, nfeatures, nclusters, SIMULATION.FILE.DIR)
  cmd.return = system(cmd, intern = F)
  stopifnot(cmd.return == 0)
}

generate.noisy.data <- function(mean, sd, outfile){
  base_data <- read.table(get.base.data.path())
  
  data <- base_data +  matrix(rnorm(n=nsamples*nfeatures, mean = mean, sd = sd), ncol=nsamples, nrow=nfeatures)
  write.table(data, file=outfile, sep = "\t", row.names = T, col.names = T)
  print(paste("Created:", outfile))
}

run.simulation <- function(files, sampling){
  dir_name <- get.clustering.results.dir.path()
  if (!dir.exists(dir_name)){
    dir.create(dir_name)
  } 
  results <- run.sampling(fnames=files[[1]], name=sampling)
  return(results)
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


generate.missing.data <- function(fraction, layer1, layer2, outdir){
  stopifnot(file.exists(layer1) & file.exists(layer2))
  data1 <- read.table(layer1)
  data2 <- read.table(layer2)
  stopifnot(all(colnames(data1) %in% colnames(data2)))
  if (! dir.exists(outdir)){
    dir.create(outdir)
  }
  # remove samples from first layer
  samples <- sample(colnames(data1), size = round((1 -fraction)*dim(data1)[2]))
  subsampled1 <- data1[,colnames(data1) %in% samples] 
  fname1 <- file.path(outdir, paste0("layer1_", fraction, ".tsv"))
  write.table(subsampled1, fname1)
  
  return(list(list(layer1=fname1, layer2=layer2)))
}

run.sampling <- function(fnames, name) {
  result_files <- list()
  subtype = name
  subtype.raw.data = get.raw.data(fnames)
  for (algorithm.name in ALGORITHM.NAMES) {
    print(paste('data', subtype, 'running algorithm', algorithm.name))
    result_dir <- get.clustering.results.dir.path()
    clustering.path = file.path(result_dir, paste0(subtype, "_", algorithm.name, ".tsv"))
    timing.path = file.path(result_dir, paste0(subtype, "_",algorithm.name, '_timing'))
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
        stop('Unknown sample order!')
      }
      result_files <- append(result_files, list(list(fname=clustering.path, subtype=subtype, algorithm=algorithm.name)))
    }
  }
  return(result_files)
}

run.evaluation <- function(results, outfile){
  algorithms <- c()
  fractions <- c()
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
      fractions <- c(fractions, r$subtype)
      vals <- c(vals, val)
      metrics <- c(metrics, metric)
    }
  }
  data <- data.frame(tool=algorithms, fraction=fractions, metric = metrics, val = vals)
  write.table(data, file=outfile, sep = "\t", row.names = F, col.names = T)
}
