job_file <- commandArgs(trailingOnly=TRUE)[1]
options(mc.cores=8)
files <- readLines("job-files.txt")

clusterMap(cl=NULL, fun=function(x) system(paste0("run.sh ", x), wait=TRUE),
  x = files, RECYCLE=TRUE, .scheduling='dynamic')


