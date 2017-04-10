library(magrittr); library(digest);
source("../../../utilities.R")

args <- commandArgs(trailingOnly=TRUE)
model_path <- args[1]
data_path <- args[2]
sub_command <- args[3]
sub_file <- args[4]

data_files <- dir(path=data_path, pattern='data-target.*rdump')

fit_targets <- list(
  `gamma-p1` = c('N', 'alpha', 'beta'),
  `weibull-p1` = c('N', 'alpha', 'beta'),
  `generalized-gamma-p1` = c('N', 'alpha', 'beta', 'nu'),
  `generalized-gamma-p2` = c('N', 'k', 'mu', 'sigma')
)

fit_target_names <- names(fit_targets)

shared_sampling_opts <- paste("method=sample",
  "algorithm=hmc engine=nuts",
  "num_samples=800 num_warmup=400 save_warmup=1 thin=1")

shared_optimizing_opts <- paste("method=optimize",
  "algorithm=lbfgs", 
  "iter=500 save_iterations=1")

shared_data_opts <- paste0("data file=")
shared_init_opts <- NULL

shared_output_opts <- paste("output refresh=10 file=")
shared_diagnostic_opts <- paste("diagnostic_file=")

stan_sample_call <- function(target, id=1, data_path, output_path, diagnostic_path) {
  o <- paste(target, shared_sampling_opts, paste0("id=", id),
    paste0(shared_data_opts, data_path),
    paste0(shared_output_opts, output_path),
    paste0(shared_diagnostic_opts, diagnostic_path))
  return(o)
}

stan_line <- list()
for ( fit in fit_target_names ) {
  for ( data in data_files ) {
    for ( id in 1:5 ) {
      data_name <- gsub(pattern='\\.rdump', x=basename(data), replacement='')
      fit_name <- paste0('model-target-', fit, '--chain-', id, '--', data_name)
      fit_hash <- digest(fit_name)
      stan_path <- file.path(model_path, fit) 
      data_full_path <- file.path(data_path, data) 
      output_full_path <- file.path('.', paste0(fit_name, '.csv')) 
      diagnostic_full_path <- file.path('.', paste0(fit_name, '--diagnostic.csv')) 
      stan_line[[fit_hash]] <- stan_sample_call(stan_path, id=id,
        data_path=data_full_path, output_path=output_full_path,
        diagnostic_path=diagnostic_full_path)
    }
  }
}

for (l in stan_line) {
  cat(paste(sub_command, l), sep="\n", file=sub_file, append=TRUE)
}

