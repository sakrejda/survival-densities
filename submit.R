library(magrittr)

job_files <-  dir(pattern='.*-job-tag.pbs')
cmd_stub_a <- "bsub -R rusage[mem=4096] -q condo_uma_nicholas_reich -W 30 -oo "
cmd_stub_b <- " -eo "
job_files_randomized <- sample(x=job_files, size=length(job_files), replace=FALSE)
job_output <- function(f) gsub(pattern=".*-chain-..-", replacement="", x=f) %>% 
  gsub(pattern="-job-tag.pbs", replacement="-output.txt")
job_error <- function(f) gsub(pattern=".*-chain-..-", replacement="", x=f) %>% 
  gsub(pattern="-job-tag.pbs", replacement="-error.txt")
for ( i in 1:length(job_files_randomized)) { 
  job <- job_files_randomized[i]; 
  cmd <- paste(cmd_stub_a, job_output(job), cmd_stub_b, job_error(job), " < ", job); 
  print(cmd) 
  system(cmd)
}

