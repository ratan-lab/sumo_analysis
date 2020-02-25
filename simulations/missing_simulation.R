source('subsampling_simulation.R')

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
