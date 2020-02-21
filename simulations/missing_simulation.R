source('subsampling_simulation.R')

generate.missing.data <- function(fraction, layer1, layer2){
  stopifnot(file.exists(layer1) & file.exists(layer2))
  data1 <- read.table(layer1)
  data2 <- read.table(layer2)
  stopifnot(all(colnames(data1) %in% colnames(data2)))
  
  # remove samples from first layer
  set.seed(RANDOM.SEED)
  samples <- sample(colnames(data1), size = round((1 -fraction)*dim(data1)[2]))
  subsampled1 <- data1[,colnames(data1) %in% samples] 
  fname1 <- file.path(SIMULATION.FILE.DIR, paste0("layer1_", fraction, ".tsv"))
  write.table(subsampled1, fname1)
  
  return(list(list(layer1=fname1, layer2=layer2)))
}
