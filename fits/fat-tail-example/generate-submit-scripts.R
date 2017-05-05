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

delays <- readRDS("delay-data.rds")
N_clean <- length(delays[['clean']])
N_contaminated <- length(delays[['contaminated']])

n_chains <- 10
chains <- 1:n_chains

data <- list(
  gamma = list(
    n_obs=N_clean,
    t_start=rep(0, N_clean), t_stop=delays[['clean']], t_truncate=rep(Inf, N_clean)
  ),
  gesgm = list(
    n_obs=N_contaminated,
    t_start=rep(0, N_contaminated), t_stop=delays[['contaminated']], t_truncate=rep(Inf, N_contaminated)
  )
)
data_names <- names(data)

model_names <- c(
  "gamma-p1",
  "generalized-gamma-p1",
  "generalized-gamma-p2",
  "gamma-exp-sum-gamma-mix-p1",
  "gamma-exp-sum-gamma-mix-p2"
)

model_data <- list(
  gamma = list(
    kappa=rep(1,2),  n_parameters=2,
    sigma_log=1/2, sigma_identity=1, sigma_logit=2/3
  ),
  generalized_gamma = list(
    kappa=rep(1,3), n_parameters=3, 
    sigma_log=1/2, sigma_identity=1, sigma_logit=2/3
  ),
  gesgm = list(
    kappa=c(10^-1, 10^-1, 1, .2, 1, 1, 1, 1), n_parameters=8, 
    sigma_log=1/2, sigma_identity=1, sigma_logit=2/3
  )
)

models_dir <- file.path("..", "..", "models",
  "stan-lang", "full")
models_dir_fix <- file.path("..", "..", "models",
  "stan-lang", "full-fix")
binary_dir <- c(models_dir, models_dir_fix)

run_data <- expand.grid(chain=chains, model=model_names, binary_dir=binary_dir, data=data_names) 

run_data <- run_data %>% mutate(
  job_id = 1:n(),
  job_label = sapply(job_id, pad, n=nchar(n())),
  chain_label = sapply(chain, pad, n=nchar(n_chains)),
  model_type = gsub(pattern='-p[0-9]', x=model, replacement=''),
  binary_type = ifelse(grepl(pattern='full-fix', x=binary_dir), 'fix', 'dev'),
  model_binary = paste0(binary_dir, model),
  model_progam = paste0(binary_dir, model, '.stan'),
  run_name = paste0(model, '_on_', data, '_with_', binary_type, '_chain_', chain_label, '_job_', job_label) %>% 
    gsub(pattern='_', x=., replacement='-'),
  data_file = paste0(model, '-on-', data, "-data.rdump"),
  init_file = paste0(model, '-on-', data, "-init.rdump"),
  output_file = paste0(run_name, "-output.csv"),
  diagnostic_file = paste0(run_name, "-diagnostic.csv")
) 

saveRDS(run_data, file='run-data.rds')


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
  data_data <- c(data[[ run_data[i, 'data'] ]], model_data[[ run_data[i, 'model_type']]])
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
    R=paste("usage[mem=4096]")
  )

  write_cmdstan_pbs(binary, args, pbs_args, prefix=bsub_prefix, pbs_file)

  my_rdump(data_data, file=run_data[i, 'data_file'])
}






