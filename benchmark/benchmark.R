# overwrite some functions
get.mkl.binary.path = function() {
  return(MKL.BINARY.PATH)
}

get.mkl.arguments.path = function() {
  return(MKL.ARGS.PATH)
}

get.mcr.root.path = function(){
  return(MCR.ROOT)
}

get.dataset.dir.path <- function() {
  return(DATASETS.PATH)
}

get.clustering.results.dir.path <- function() {
  if (!dir.exists(RESULTS.DIR.PATH)){
    dir.create(RESULTS.DIR.PATH)
  }
  return(RESULTS.DIR.PATH)
}

get.clinical.params.dir <- function() {
  return(CLINICAL.PARAMS.DIR)
}

get.tables.dir.path <- function() {
  results.dir.path = get.clustering.results.dir.path()
  tables.dir = file.path(results.dir.path, 'tables')
  if (!dir.exists(tables.dir)){
    dir.create(tables.dir)
  }
  return(tables.dir)
}

set.omics.list.attr <- function(subtype.raw.data, subtype.data) {
  attr(subtype.raw.data[[1]], 'is.seq') = subtype.data$is.rna.seq
  attr(subtype.raw.data[[1]], 'is.met') = F
  attr(subtype.raw.data[[2]], 'is.seq') = F
  attr(subtype.raw.data[[2]], 'is.met') = T
  attr(subtype.raw.data[[3]], 'is.seq') = subtype.data$is.mirna.seq
  attr(subtype.raw.data[[3]], 'is.met') = F
  return(subtype.raw.data)
}

get.vars.path <- function(){
  return(VARS.FNAME)
}

load.libraries <- function() {
  #tools
  library('PMA')
  library('R.matlab')
  library('SNFtool')
  library('PINSPlus')
  library('LRAcluster')
  reticulate::use_python(Sys.which('python3'), required = TRUE)
  library('reticulate')
  library("DESeq2")
  stopifnot(file.exists(get.vars.path()))
  system(paste("sh", get.vars.path()))
  #analysis
  library('survival')
}

run.mkl <- function(omics.list, subtype.data) {
  start = Sys.time()
  omics.list = log.and.normalize(omics.list, subtype.data)
  subtype = subtype.data$name
  omics.list = lapply(omics.list, normalize.matrix)
  time.taken = as.numeric(Sys.time() - start, units='secs')
  export.subtype.to.mkl(omics.list, subtype)
  
  start = Sys.time()
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
  time.taken2 = as.numeric(Sys.time() - start, units='secs')
  clustering = get.mkl.clustering(subtype)
  return(list(clustering=clustering, 
              timing=time.taken + time.taken2))
}

get.mkl.clustering <- function(dir.name) {
  folder.path = file.path(get.mkl.arguments.path(), dir.name)
  output.path = file.path(folder.path, 'output')
  output.files = list.files(output.path)
  clustering = read.csv(file.path(output.path, output.files[grep('clusters', output.files)]))[,2]
  return(clustering)
}

export.subtype.to.mkl <- function(omics.list, dir.name, run.id=NULL) {
  folder.path = file.path(get.mkl.arguments.path(), dir.name, toString(run.id))

  if (!dir.exists(folder.path)) {
    dir.create(folder.path)
  }

  kernels.path = file.path(folder.path, 'kernels')

  if (!dir.exists(kernels.path)) {
    dir.create(kernels.path)
  }
  clear.dir(kernels.path)

  gammas = 10 ** seq(-6, 6, by=3)
  for (i in 1:length(omics.list)) {
    for (j in 1:length(gammas)) {
      datum = omics.list[[i]]
      gamma = gammas[[j]] / (2*nrow(datum)**2)
      mat = radial.basis(datum, gamma)
      R.matlab::writeMat(file.path(kernels.path, paste(i, '_', j, sep='')), mat=mat)
    }
  }

  output.path = file.path(folder.path, 'output')
  if (!dir.exists(output.path)) {
    dir.create(output.path)
  }
  clear.dir(output.path)

  write.table(colnames(omics.list[[1]]), file=file.path(folder.path, 'ids'),
              quote=F, row.names = F, col.names = F)

}

benchmark.omics.num.clusters <- function(benchmark.results, omics='all') {
  num.clusters = matrix(1, ncol=length(SUBTYPES.DATA), nrow=length(ALGORITHM.NAMES))
  rownames(num.clusters) = ALGORITHM.NAMES
  colnames(num.clusters) = sapply(SUBTYPES.DATA, function(x) x$name)
  all.clusterings = benchmark.results$all.clusterings
  for (i in 1:length(all.clusterings)) {
    subtype = colnames(num.clusters)[i]
    subtype.clusterings = all.clusterings[[subtype]]
    for (j in 1:length(subtype.clusterings)) {
      clustering = subtype.clusterings[[ALGORITHM.NAMES[j]]][[omics]]
      num.clusters[j, i] = as.numeric(max(clustering))
    }
  }
  return(num.clusters)
}

get.clinical.params <- function(subtype.name) {
  clinical.data.path = file.path(get.clinical.params.dir(), subtype.name)
  clinical.params = read.table(clinical.data.path, sep='\t', header=T, row.names = 1, stringsAsFactors = F)
  rownames.with.duplicates = get.fixed.names(rownames(clinical.params))  
  clinical.params = clinical.params[!duplicated(rownames.with.duplicates),]
  rownames(clinical.params) = rownames.with.duplicates[!duplicated(rownames.with.duplicates)]
  return(clinical.params)
}

analyze.benchmark <- function() {
  all.clusterings = list()
  all.timings = list()
  for (i in 1:length(SUBTYPES.DATA)) {
    current.subtype.data = SUBTYPES.DATA[[i]]
    subtype = current.subtype.data$name
    subtype.raw.data = get.raw.data(subtype, 
                                    only.primary=current.subtype.data$only.primary)
    
    all.clusterings[[subtype]] = list()
    all.timings[[subtype]] = list()
    
    for (algorithm.name in ALGORITHM.NAMES) {
      all.clusterings[[subtype]][[algorithm.name]] = list()
      all.timings[[subtype]][[algorithm.name]] = list()
      for (j in names(OMIC.SUBSETS)) {
        clustering.path = file.path(get.clustering.results.dir.path(),
                                    paste(subtype, algorithm.name, j, sep='_'))
        timing.path = file.path(get.clustering.results.dir.path(),
                                paste(subtype, algorithm.name, j, 'timing', sep='_'))
        load(clustering.path)
        load(timing.path)
        if (!any(is.na(clustering))) {
          names(clustering) = colnames(subtype.raw.data[[1]])
        }
        
        all.clusterings[[subtype]][[algorithm.name]][[j]] = clustering
        all.timings[[subtype]][[algorithm.name]][[j]] = timing
      }
    }
  }
  return(list(all.clusterings=all.clusterings, all.timings=all.timings))
}

run.benchmark <- function() {
  for (i in 1:length(SUBTYPES.DATA)) {
    current.subtype.data = SUBTYPES.DATA[[i]]
    subtype = current.subtype.data$name
    subtype.raw.data = get.raw.data(subtype,
                                    only.primary=current.subtype.data$only.primary)

    subtype.raw.data = set.omics.list.attr(subtype.raw.data,
                                           current.subtype.data)

    for (algorithm.name in ALGORITHM.NAMES) {
      for (j in names(OMIC.SUBSETS)) {
        set.seed(42)
        print(paste('subtype', subtype, 'running algorithm', algorithm.name, j))
        clustering.path = file.path(get.clustering.results.dir.path(),
                                    paste(subtype, algorithm.name, j, sep='_'))
        timing.path = file.path(get.clustering.results.dir.path(),
                                paste(subtype, algorithm.name, j, 'timing', sep='_'))

        if (!file.exists(clustering.path)) {
          algorithm.func.name = paste0('run.', algorithm.name)
          algorithm.func = get(algorithm.func.name)
          if (j == 'all') {
            cur.iteration.data = subtype.raw.data
          } else {
            cur.iteration.data = subtype.raw.data[as.numeric(j)]
          }
          algorithm.ret = algorithm.func(cur.iteration.data, current.subtype.data)
          clustering = algorithm.ret$clustering
          timing = algorithm.ret$timing
          print('before saving')
          save(clustering, file = clustering.path)
          save(timing, file = timing.path)
        }
      }
    }
  }
}

get.empirical.clinical <- function(clustering, clinical.values, is.chisq) {
  set.seed(42)
  clustering <- as.factor(clustering)
  if (is.chisq) {
    clustering.with.clinical = cbind(clustering, clinical.values)
    tbl = table(as.data.frame(clustering.with.clinical))
    test.res = chisq.test(tbl)
  } else {
    test.res = kruskal.test(as.numeric(clinical.values), clustering)
  }
  orig.pvalue = test.res$p.value
  num.iter = 1000
  total.num.iters = 0
  total.num.extreme = 0
  should.continue = T

  while (should.continue) {
    print('another iteration in empirical clinical')
    perm.pvalues = as.numeric(mclapply(1:num.iter, function(i) {
      cur.clustering = sample(clustering)
      names(cur.clustering) = names(clustering)

      if (is.chisq) {
        clustering.with.clinical = cbind(cur.clustering, clinical.values)
        tbl = table(as.data.frame(clustering.with.clinical))
        test.res = chisq.test(tbl)
      } else {
        test.res = kruskal.test(as.numeric(clinical.values), cur.clustering)
      }
      cur.pvalue = test.res$p.value
      return(cur.pvalue)
    }, mc.cores=get.mc.cores()))
    total.num.iters = total.num.iters + num.iter
    total.num.extreme = total.num.extreme + sum(perm.pvalues <= orig.pvalue)

    binom.ret = binom.test(total.num.extreme, total.num.iters)
    cur.pvalue = binom.ret$estimate
    cur.conf.int = binom.ret$conf.int

    sig.threshold = 0.05
    is.threshold.in.conf = cur.conf.int[1] < sig.threshold & cur.conf.int[2] > sig.threshold
    if (!is.threshold.in.conf | total.num.iters > 1e5) {
      should.continue = F
    }
  }
  return(cur.pvalue)
}

get.empirical.surv <- function(clustering, subtype) {
  set.seed(42)
  surv.ret = check.survival(clustering, subtype)
  orig.chisq = surv.ret$chisq
  orig.pvalue = get.logrank.pvalue(surv.ret)
  # The initial number of permutations to run
  num.perms = round(min(max(10 / orig.pvalue, 1000), 1e6))
  should.continue = T

  total.num.perms = 0
  total.num.extreme.chisq = 0

  while (should.continue) {
    print('Another iteration in empirical survival calculation')
    print(num.perms)
    perm.chisq = as.numeric(mclapply(1:num.perms, function(i) {
      cur.clustering = sample(clustering)
      names(cur.clustering) = names(clustering)
      cur.chisq = check.survival(cur.clustering, subtype)$chisq
      return(cur.chisq)
    }, mc.cores=get.mc.cores()))

    total.num.perms = total.num.perms + num.perms
    total.num.extreme.chisq = total.num.extreme.chisq + sum(perm.chisq >= orig.chisq)

    binom.ret = binom.test(total.num.extreme.chisq, total.num.perms)
    cur.pvalue = binom.ret$estimate
    cur.conf.int = binom.ret$conf.int

    print(c(total.num.extreme.chisq, total.num.perms))
    print(cur.pvalue)
    print(cur.conf.int)

    sig.threshold = 0.05
    is.conf.small = ((cur.conf.int[2] - cur.pvalue) < min(cur.pvalue / 10, 0.01)) & ((cur.pvalue - cur.conf.int[1]) < min(cur.pvalue / 10, 0.01))
    is.threshold.in.conf = cur.conf.int[1] < sig.threshold & cur.conf.int[2] > sig.threshold
    if ((is.conf.small & !is.threshold.in.conf) | (total.num.perms > 2e7)) {
      should.continue = F
    } else {
      num.perms = 1e5
    }
  }

  return(list(pvalue = cur.pvalue, conf.int = cur.conf.int, total.num.perms=total.num.perms,
              total.num.extreme.chisq=total.num.extreme.chisq))
}

get.mc.cores <- function(){
  return(MC.CORES)
}

perform.all.analyses <- function(benchmark.ret) {
  
  for (i in 1:4) {
    
    cur.func = list(benchmark.omics.time, benchmark.omics.num.clusters,
                    benchmark.omics.surv, benchmark.omics.clinical)[[i]]
    for (omic.subset in names(OMIC.SUBSETS)) {
      print(paste('current omic subset ', omic.subset))
      benchmark.data = cur.func(benchmark.ret, omic.subset)
      
      displayed.benchmark.data = benchmark.data
      colnames(displayed.benchmark.data)[1:ncol(displayed.benchmark.data)] = 
        sapply(as.list(colnames(displayed.benchmark.data)[1:ncol(displayed.benchmark.data)]), 
               subtype.to.display.name)
      rownames(displayed.benchmark.data) = unlist(ALGORITHM.DISPLAY.NAMES[rownames(displayed.benchmark.data)])
      print.matrix.latex.format(displayed.benchmark.data)
      
      table.name = c('runtime', 'num_cluster', 'survival', 'clinical')[i]
      write.csv(displayed.benchmark.data, file=file.path(get.tables.dir.path(), paste0(table.name, '_', OMIC.SUBSETS[[omic.subset]], '.csv')))
    }
    
    
    print('------------------------')
  }
}



plot_results <- function(){

  library(tidyverse)

  clin <- read_csv(file.path(get.tables.dir.path(), "clinical_multi_omics.csv"))
  surv <- read_csv(file.path(get.tables.dir.path(), "survival_multi_omics.csv"))

  clin_long <- clin %>%
    gather(cancer, clin, -X1) %>%
    mutate(tool=X1) %>%
    select(-X1)

  surv %>%
    gather(cancer, surv, -X1) %>%
    mutate(tool=X1) %>%
    select(-X1) %>%
    full_join(clin_long) %>%
    ggplot() + geom_point(aes(x=surv, y=clin, color=tool)) + facet_wrap(cancer~., scales="free") +
    xlab("-log10(logrank pvalue)") + ylab("# enriched clinical parameters") + labs(color="") +
    theme(legend.position="bottom") + geom_vline(xintercept = -log10(0.05), color="red")

  p <- last_plot()
  return(p)

}
