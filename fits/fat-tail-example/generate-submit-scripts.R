library(stannis)

my_rdump <- function(x, file) {
  rstan::stan_rdump(names(x), file, envir=list2env(x))
}

library(rstan); library(dplyr); library(magrittr)

delays <- readRDS("delay-data.rds")
N_clean <- length(delays[['clean']])
N_contaminated <- length(delays[['contaminated']])

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

model_inits <- list(
  gamma_on_gamma = NULL,
  gamma_on_gesgm = NULL,
  generalized_gamma_on_gamma = NULL,
  generalized_gamma_on_gesgm = NULL,
  gesgm_on_gamma = list(
    mu_0 = sqrt(3), sigma_0=sqrt(sqrt(3)), q_0=-5, delta_0=0
  ),
  gesgm_on_gesgm = list(
    mu_0 =log(sqrt(3)*sqrt(3)), sigma_0=log(sqrt(3)), q_0=-1, delta_0=log(5/1)-3
  )
)

data_names <- names(data)
model_data_names <- names(model_data) 

run_data <- expand.grid(model=model_data_names, data=data_names) 
run_names <- run_data %>% as.matrix %>% apply(1, paste, collapse='_on_') %>%
  gsub(pattern='_', x=., replacement='-')
run_data <- cbind(run_data, run_names) %>% data.frame(check.names=FALSE) %>%
  mutate(data_files = paste0(run_names, "-data.rdump")) %>%
  mutate(init_files = paste0(run_names, "-init.rdump")) %>%
  mutate(output_files = paste0(run_names, "-X-output.csv")) %>%
  mutate(diagnostic_files = paste0(run_names, "-X-diagnostic.csv")) 

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

models_dir <- file.path("..", "..", "models",
  "stan-lang", "full")

models <- c(
  gamma = "gamma-p1",
  generalized_gamma = "generalized-gamma-p1",
  gesgm = "gamma-exp-sum-gamma-mix-p2"
)

model_binaries <- file.path(models_dir, models) %>% `names<-`(names(models))
model_programs <- paste(model_binaries, 'stan', sep='.')

saveRDS(model_binaries, 'model-binary-files.rds')
saveRDS(model_programs, 'model-program-files.rds')


n_chains <- 10
n_run_types <- nrow(run_data)
n_runs <- n_run_types*n_chains

pad <- function(x, n) {
  nc <- nchar(as.character(x))
  if ((n-nc) < 0)
    n <- nc
  pads <- rep('0', n-nc) %>% paste(collapse='')
  if ((n-nc) > 0) 
    x <- paste0(pads, x, collapse='', sep='')
  return(x)
}

width <- nchar(as.character(n_runs))
for (r in 1:n_runs) {
  chain_id <- r
  chain_label <- pad(chain_id, width)
  i <- ((r-1) %% n_run_types) + 1
  binary <- model_binaries[[ run_data[i,'model'] ]] 
  data_data <- c(data[[ run_data[i, 'data'] ]], model_data[[ run_data[i, 'model']]])
  files <- list(
    data = run_data[i, 'data_files'],
    init = run_data[i, 'init_files'],
    output = gsub(pattern='-X-', x=run_data[i, 'output_files'], replacement=paste0('-',chain_label, '-')),
    diagnostic = gsub(pattern='-X-', x=run_data[i, 'diagnostic_files'], replacement=paste0('-',chain_label, '-'))
  )
  if (is.null(model_inits[[ run_data[i, 'run_names'] ]])) 
    files[['init']] <- NULL
  binary <- paste0("  cd $PBS_O_WORKDIR \n", model_binaries[ run_data[i, 'model'] ] )
  args <- stan_args(chain_label, files)
  pbs_file <- paste0(run_data[i, 'run_names'], '-id-', chain_label, '-job.pbs')
  pbs_args <- list(
    N=paste0(run_data[i, 'run_names'], '-id-', chain_label, '-job-tag')
  )
  write_cmdstan_pbs(binary, args, pbs_args, pbs_prefix, pbs_file)

  my_rdump(data_data, file=run_data[i, 'data_files'])
  if (!is.null(model_inits[[ run_data[i, 'run_names'] ]] ))
    my_rdump(model_inits[[ run_data[i, 'run_names'] ]], file=run_data[i, 'init_files'])
}






