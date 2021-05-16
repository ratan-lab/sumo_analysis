# overwrite some functions
get.random.seed <- function(){
  return(RANDOM.SEED)
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
  library('iClusterPlus')
  library("CIMLR")
  reticulate::use_python(Sys.which('python3'), required = TRUE)
  library('reticulate')
  library("DESeq2")
  library('cluster')
  stopifnot(file.exists(get.vars.path()))
  system(paste("sh", get.vars.path()))
  #analysis
  library('survival')
}

run.iCluster <- function(omics.list, subtype.data) {
  omics.list = log.and.normalize(omics.list, subtype.data, normalize = F)
  
  start = Sys.time()
  subtype = subtype.data$name
  dev.ratios = c()
  icluster.rets = list()
  
  if (length(omics.list) == 1) {
    icluster.ret = iClusterPlus::tune.iClusterBayes(cpus=get.mc.cores(), t(omics.list[[1]]), 
                                                    K=1:(MAX.NUM.CLUSTERS - 1), type=c('gaussian'))$fit
  } else {
    icluster.ret = iClusterPlus::tune.iClusterBayes(cpus=get.mc.cores(), t(omics.list[[1]]), 
                                                    t(omics.list[[2]]), 
                                                    t(omics.list[[3]]), 
                                                    K=1:(MAX.NUM.CLUSTERS - 1), type=rep('gaussian', 3))$fit
  }
  dev.ratios = lapply(1:(MAX.NUM.CLUSTERS - 1), function(i) icluster.ret[[i]]$dev.ratio)
  
  print('dev.ratios are:')
  print(dev.ratios)
  
  optimal.solution = icluster.ret[[which.max(dev.ratios)]]
  time.taken = as.numeric(Sys.time() - start, units='secs')
  return(list(clustering=optimal.solution$clusters, 
              timing=time.taken))
}

run.nemo <- function(omics.list, subtype.data, num.clusters=NULL, is.missing.data=F, num.neighbors=NA) {
  omics.list = log.and.normalize(omics.list, subtype.data)
  start = Sys.time()
  clustering = nemo.clustering(omics.list, num.clusters, num.neighbors)
  time.taken = as.numeric(Sys.time() - start, units='secs')
  return(list(clustering=clustering, timing=time.taken))
}

run.mcca <- function(omics.list, subtype.data, known.num.clusters=NULL, penalty=NULL, rep.omic=1) {
  
  if (length(omics.list) == 1) {
    return(list(clustering=rep(1, ncol(omics.list[[1]])), timing=1))
  }
  omics.list = log.and.normalize(omics.list, subtype.data, 
                                 normalize = T,
                                 filter.var = T)
  start = Sys.time()
  max.dim = min(MAX.NUM.CLUSTERS, nrow(omics.list[[1]]))
  
  subtype = subtype.data$name
  omics.transposed = lapply(omics.list, t)
  
  cca.ret = PMA::MultiCCA(omics.transposed, 
                          ncomponents = max.dim, penalty=penalty)
  sample.rep = omics.transposed[[rep.omic]] %*% cca.ret$ws[[rep.omic]]
  
  explained.vars = sapply(1:max.dim, 
                          function(i) sum(unlist(apply(sample.rep[1:i,,drop=F], 2, var))))
  
  dimension = get.elbow(explained.vars, is.max=F)
  print(dimension)
  sample.rep = sample.rep[,1:dimension]
  sils = c()
  clustering.per.num.clusters = list()
  if (is.null(known.num.clusters)) {
    for (num.clusters in 2:max.dim) {
      cur.clustering = kmeans(sample.rep, num.clusters, iter.max=100, nstart=30)$cluster  
      sil = get.clustering.silhouette(list(t(sample.rep)), cur.clustering)
      sils = c(sils, sil)
      clustering.per.num.clusters[[num.clusters - 1]] = cur.clustering
    }
    cca.clustering = clustering.per.num.clusters[[which.min(sils)]]
  } else {
    cca.clustering = kmeans(sample.rep, known.num.clusters, iter.max=100, nstart=30)$cluster  
  }
  
  time.taken = as.numeric(Sys.time() - start, units='secs')
  return(list(clustering=cca.clustering, timing=time.taken))
}

run.snf <- function(omics.list, subtype.data) {
  omics.list = log.and.normalize(omics.list, subtype.data)
  start = Sys.time()
  subtype = subtype.data$name
  alpha=0.5
  T.val=30
  num.neighbors = round(ncol(omics.list[[1]]) / 10)
  similarity.data = lapply(omics.list, function(x) {affinityMatrix(dist2(as.matrix(t(x)),as.matrix(t(x))), 
                                                                   num.neighbors, alpha)})
  if (length(similarity.data) == 1) {
    W = similarity.data[[1]]
  } else {
    W = SNF(similarity.data, num.neighbors, T.val)  
  }
  
  num.clusters = estimateNumberOfClustersGivenGraph(W, 2:MAX.NUM.CLUSTERS)[[3]]  
  clustering = spectralClustering(W, num.clusters)
  time.taken = as.numeric(Sys.time() - start, units='secs')
  return(list(clustering=clustering, timing=time.taken))
}

run.pins <- function(omics.list, subtype.data) {
  omics.list = log.and.normalize(omics.list, subtype.data, normalize = F)
  start = Sys.time()
  subtype = subtype.data$name
  omics.transposed = lapply(omics.list, t)
  if (length(omics.list) == 1) {
    pins.ret = PINSPlus::PerturbationClustering(data=omics.transposed[[1]], kMax = MAX.NUM.CLUSTERS,
                                                ncore = get.mc.cores())
    clustering = pins.ret$cluster
    
  } else {
    pins.ret = PINSPlus::SubtypingOmicsData(dataList=omics.transposed, kMax = MAX.NUM.CLUSTERS,
                                            ncore = get.mc.cores())
    clustering = pins.ret$cluster2
  }
  time.taken = as.numeric(Sys.time() - start, units='secs')
  return(list(clustering=clustering, timing=time.taken))
}

run.cimlr <- function(omics.list, subtype.data) {
  cimlr_preprocess <- function(omic, is.seq, max_threshold=10, min_threshold=-10){
    if (is.seq){
      omic = log(1+omic)
      omic = normalize.matrix(omic)
      # z-score filtering
      omic[omic > max_threshold] = max_threshold
      omic[omic < min_threshold] = min_threshold
      # normalize data to the [0,1] range
      max_threshold <- max(omic)
      min_threshold <- min(omic)
      omic <- (omic - min_threshold + .Machine$double.eps)/(max_threshold - min_threshold + .Machine$double.eps)
    }
    return(omic)
  }
  
  omics.list = lapply(omics.list, function(x){cimlr_preprocess(omic=x, is.seq=attr(x, 'is.seq'))})
  start = Sys.time()

  subtype = subtype.data$name
  NUMC = 2:MAX.NUM.CLUSTERS
  cores.ratio = ifelse(MC.CORES / detectCores() <= 1, MC.CORES / detectCores(), 1)
  estimate_k <- CIMLR_Estimate_Number_of_Clusters(omics.list, NUMC = NUMC, cores.ratio = cores.ratio)
  print(paste0("Best number of clusters, K1 heuristic: ",NUMC[which.min(estimate_k$K1)],", K2 heuristic: ",NUMC[which.min(estimate_k$K2)]))
  num.clusters <- NUMC[which.min(estimate_k$K1)]
  
  res <- CIMLR(omics.list, c=num.clusters, cores.ratio = cores.ratio)
  clustering <- res$y$cluster

  time.taken = as.numeric(Sys.time() - start, units='secs')
  return(list(clustering=clustering, timing=time.taken))
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
        set.seed(get.random.seed())
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
  set.seed(get.random.seed())
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
  set.seed(get.random.seed())
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
