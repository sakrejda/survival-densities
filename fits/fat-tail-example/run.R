library(parallel)
cl <- makeForkCluster(nnodes=8)
job_file <- commandArgs(trailingOnly=TRUE)[1]
source_dir <- commandArgs(trailingOnly=TRUE)[2]
files <- readLines("job-files.txt")

clusterMap(cl=cl, fun=function(x, source_dir) system(paste0(file.path(source_dir,"run.sh"), " ", x), wait=TRUE),
  x = files, source_dir=source_dir, RECYCLE=TRUE, .scheduling='dynamic')


