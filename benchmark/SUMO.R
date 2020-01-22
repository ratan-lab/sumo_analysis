library(reticulate)
# NOTE when using reticulate you may encounter an error in initialize_python which prevents python bindings being loaded,
# to solve this issue make sure "--enable-shared" option was used when building from source your default python version available in sumo_analysis directory

sumo.preprocess.data <- function(data){
  is.seq <- attr(data, "is.seq")
  is.met <- attr(data, "is.met")
  if (is.seq){ # not methylation and not microarray data, which is previously normalized
    # normalization 
    dds <- DESeqDataSetFromMatrix(round(data), colData = data.frame(sample = colnames(data)), design = ~1)
    vsd <- varianceStabilizingTransformation(dds, blind = TRUE)
    data <- assay(vsd)
  } else {
    if (is.met){ #methylation
      stopifnot(max(data)<= 1 & min(data)>=0) #beta values
      beta <- as.matrix(data)
      Mvals <- log2(beta/(1-beta))
      data <- Mvals
    } 
  }
  # standarization
  standard.data <- t(scale(t(as.matrix(data))))
  print(dim(standard.data))
  return(standard.data)
}

get.sumo.files.dir <- function(){
  return(SUMO.FILES.DIR)
}

get.sumo.path <- function(){
  return(SUMO.PATH)
}

sumo.clustering <- function(k, fname, outdir, mc.cores = get.mc.cores()){
  log <- paste0(paste(strsplit(fname, "\\.")[[1]][1:2], collapse = '.'), ".log")
  command <- paste(get.sumo.path(), "run","-t", mc.cores,"-log DEBUG", "-logfile", log, fname, k, outdir)
  print(command)
  command.return = system(command)
  stopifnot(command.return == 0)

  #TODO select the best results and point to them below
  results_file <- file.path(outdir, paste0("k",2), "sumo_results.npz")
  np <- import("numpy")
  npz <- np$load(results_file, allow_pickle = T)
  clusters <- npz$f[["clusters"]]
  
  clustering <- unlist(clusters[,2])
  names(clustering) <- unlist(clusters[,1])
  return(clustering)
}

sumo.prepare <- function(omics.list, subtype.data, prepare_files, outfile){
  fnames = NULL
  for (i in prepare_files){
    if (is.null(fnames)){
      fnames <- i
    } else {
      fnames <- paste(fnames, i, sep=",")
    }
  }
  plot <- paste0(paste(strsplit(outfile, "\\.")[[1]][1:2], collapse = '.'), ".png")
  log <- paste0(paste(strsplit(outfile, "\\.")[[1]][1:2], collapse = '.'), ".log")

  command <- paste(get.sumo.path(), "prepare", "-plot", plot, "-logfile", log, "-atol", 0.1, fnames, outfile)
  print(command)
  command.return = system(command)
  stopifnot(command.return == 0)
}

run.sumo <- function(omics.list, subtype.data, num.clusters=NULL, mc.cores=get.mc.cores(), file_dir= get.sumo.files.dir()){
  start = Sys.time()
  if (!dir.exists(file_dir)){
    dir.create(file_dir)
  }
  
  prepare_files <- c()
  for (i in 1:length(omics.list)){
    fname <- file.path(get.sumo.files.dir(), paste0(subtype.data$name, "_", i, ".txt"))
    if (!file.exists(fname)){
      omic <- sumo.preprocess.data(omics.list[[i]])
      omics.list[[i]] <- omic
      write.table(omic, fname, sep="\t", na="NA")
    }
    prepare_files <- c(prepare_files, fname)
  }
  
  outfile = file.path(file_dir, paste0("prepared.", subtype.data$name, ".npz"))
  sumo.prepare(omics.list, subtype.data, prepare_files, outfile=outfile)
  stopifnot(file.exists(outfile))
  
  if(is.null(num.clusters)){
    num.clusters=paste(2,MAX.NUM.CLUSTERS, sep=",")
  }
  
  outdir = file.path(file_dir, subtype.data$name)
  clustering = sumo.clustering(num.clusters, outfile, outdir, mc.cores)
  time.taken = as.numeric(Sys.time() - start, units='secs')
  return(list(clustering=clustering, timing=time.taken))
}
