library(magrittr)

args <- commandArgs(trailingOnly=TRUE)
model_path <- args[1]
data_path <- args[2]
output_path <- args[3]

source("../utilities.R")

stan_models <- c('gamma-p1.stan', 'weibull-p1.stan',
  'generalized-gamma-p1.stan', 'generalized-gamma-p2.stan')

data_sets <- c('gamma-p1.rds', 'weibull-p1.rds',
  'generalized-gamma-p1.rds')


targets <- expand.grid(model=stan_models, data=data_sets) %>% apply(1, list)
n_targets <- length(targets)

for (i in 1:n_targets) {
  model_name <- drop_extension(targets[[i]][[1]]['model']) 
  data_name <- drop_extension(targets[[i]][[1]]['data']) 
  target_name <- paste(model_name, data_name, sep='--')
  model_full_path <- file.path(model_path, paste0(model_name,'.stan'))
  data_full_path <- file.path(data_path, paste0(data_name,'.rds'))
  output_full_path <- file.path(output_path, paste0(target_name,'.rds'))
  fit(model_full_path, data_full_path, output_full_path)
}


