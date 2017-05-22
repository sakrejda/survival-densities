library(stannis)

my_rdump <- function(x, file) {
  rstan::stan_rdump(names(x), file, envir=list2env(x))
}

library(rstan); library(dplyr); library(magrittr)

pad <- function(x, n) {
  nc <- nchar(as.character(x))
  if ((n-nc) < 0)
    n <- nc
  pads <- rep('0', n-nc) %>% paste(collapse='')
  if ((n-nc) > 0) 
    x <- paste0(pads, x, collapse='', sep='')
  return(x)
}

n_chains <- 10
chains <- 1:n_chains

### Data:
delays <- readRDS("delays.rds")
N_gamma <- lapply(delays, function(x) length(x[['gamma']]))
N_gesgm <- lapply(delays, function(x) length(x[['gesgm']]))
n_data_sets <- length(delays) 
data_sets <- 1:n_data_sets 
data_types <- 2

data <- mapply(function(delays, N_gamma, N_gesgm) list(
  gamma = list(
    n_obs=N_gamma,
    t_start=rep(0, N_gamma), t_stop=delays[['gamma']], t_truncate=rep(Inf, N_gamma)
  ),
  gesgm = list(
    n_obs=N_gesgm,
    t_start=rep(0, N_gesgm), t_stop=delays[['gesgm']], t_truncate=rep(Inf, N_gesgm)
  )
), delays=delays, N_gamma=N_gamma, N_gesgm=N_gesgm, SIMPLIFY=FALSE)
data_names <- c("gamma", "gesgm")

### Models:
model_names <- c(
  "gamma-p1",
  "generalized-gamma-p1",
  "generalized-gamma-p2",
  "gamma-exp-sum-gamma-mix-p1",
  "gamma-exp-sum-gamma-mix-p2"
)

model_data <- list(
  gamma_p1 = list(
    kappa=rep(1,2),  n_parameters=2,
    sigma_log=1/2, sigma_identity=1, sigma_logit=2/3
  ),
  generalized_gamma_p1 = list(
    kappa=rep(1,3), n_parameters=3, 
    sigma_log=1/2, sigma_identity=1, sigma_logit=2/3
  ),
  generalized_gamma_p2 = list(
    kappa=rep(1,3), n_parameters=3, 
    sigma_log=1/2, sigma_identity=1, sigma_logit=2/3
  ),
  gamma_exp_sum_gamma_mix_p1 = list(
    kappa=c(10^-1, 10^-1, 1, .2, 1, 1, 1, 1), n_parameters=8, 
    sigma_log=1/2, sigma_identity=1, sigma_logit=2/3
  ),
  gamma_exp_sum_gamma_mix_p2 = list(
    kappa=c(1, 1, 1, .2, 1, 1, 1, 1), n_parameters=8, 
    sigma_log=1/2, sigma_identity=1, sigma_logit=2/3
  )
)

# Binaries
models_dir_fix <- file.path("..", "..", "models",
  "stan-lang", "full-fix")
binary_dir <- models_dir_fix


# Grid of runs:
run_data <- expand.grid(chain=chains, model=model_names, binary_dir=binary_dir, data_type=data_names, data_sets=data_sets) 

run_data <- run_data %>% mutate(
  job_id = 1:n(),
  job_label = sapply(job_id, pad, n=nchar(n())),
  chain_label = sapply(chain, pad, n=nchar(n_chains)),
  model_type = gsub(pattern='-p[0-9]', x=model, replacement=''),
  binary_type = ifelse(grepl(pattern='full-fix', x=binary_dir), 'fix', 'dev'),
  model_binary = file.path(binary_dir, model),
  model_progam = paste0(model_binary, '.stan'),
  run_name = paste0(model, '_on_', data_type, '_with_', binary_type, '_chain_', chain_label, '_job_', job_label) %>% 
    gsub(pattern='_', x=., replacement='-'),
  data_file = paste0(model, '-on-', data_type, "-job-", job_label, "-data.rdump"),
  init_file = paste0(model, '-on-', data_type, "-init.rdump"),
  output_file = paste0(run_name, "-output.csv"),
  diagnostic_file = paste0(run_name, "-diagnostic.csv")
) 

saveRDS(run_data, file='run-data.rds')

# Construct and save run files:
stan_args <- function(chain_id, files) {
  o <- paste0(
    "method=sample algorithm=hmc engine=nuts ", 
    "num_samples=500 num_warmup=1000 save_warmup=1 thin=1 ",
    "id=", chain_id, " ",
    "data file=", files[['data']], " "
  )
  if (!is.null(files[['init']])) {
    o <- paste0(o, "init=", files[['init']], " ")
  }
  o <- paste0(o, "output refresh=10 ",
    "file=", files[['output']]," ", 
    "diagnostic_file=", files[['diagnostic']], " "
  )
  return(o)
}


saveRDS(run_data[['model_binary']], 'model-binary-files.rds')
saveRDS(run_data[['model_program']], 'model-program-files.rds')


for (i in 1:nrow(run_data)) {
  chain_label <- run_data[i, 'chain_label']
  binary <- run_data[i,'model_binary'] 
  data_data <- c(data[[run_data[i, 'data_sets'] ]][[run_data[i, 'data_type']]], 
    model_data[[ gsub(pattern='-', replacement='_', x=run_data[i, 'model'])]])
  files <- list(
    data = run_data[i, 'data_file'],
    init = NULL,  ## We don't want to require inits.
    output = run_data[i, 'output_file'],
    diagnostic = run_data[i, 'diagnostic_file'] 
  )
  args <- stan_args(chain_label, files)
  pbs_file <- paste0(run_data[i, 'run_name'], '-job-tag.pbs')
  pbs_args <- list(
    J=paste0(run_data[i, 'run_name'], '-job-tag'),
    oo=paste0(run_data[i, 'run_name'], '-terminal-output.txt'),
    eo=paste0(run_data[i, 'run_name'], '-terminal-errors.txt'),
    W=paste0("00:30"),
    q=paste0("condo_uma_nicholas_reich"),
    R=paste("rusage[mem=4096]")
  )

  write_cmdstan_pbs(binary, args, pbs_args, prefix=bsub_prefix, pbs_file)

  my_rdump(data_data, file=run_data[i, 'data_file'])
}






