my_rdump <- function(x, file) {
  rstan::stan_rdump(names(x), file, envir=list2env(x))
}

library(rstan); library(dplyr); library(magrittr)

delays <- readRDS("delay-data.rds")
N_clean <- length(delays[['clean']])
N_contaminated <- length(delays[['contaminated']])

gamma_on_gamma_data <- list(
  n_obs=N_clean,
  t_start=rep(0, N_clean), t_stop=delays[['clean']], t_truncate=rep(Inf, N_clean),
  kappa=rep(1,2),  n_parameters=2,
  sigma_log=1/2, sigma_identity=1, sigma_logit=2/3
) %>% my_rdump(file='gamma-on-gamma-data.rdump')
  
gamma_on_gesgm_data <- list(
  n_obs=N_contaminated,
  t_start=rep(0, N_contaminated), t_stop=delays[['contaminated']], t_truncate=rep(Inf, N_contaminated),
  kappa=rep(1,2), n_parameters=2,
  sigma_log=1/2, sigma_identity=1, sigma_logit=2/3
) %>% my_rdump(file='gamma-on-gesgm-data.rdump')

generalized_gamma_on_gamma_data <- list(
  n_obs=N_clean, 
  t_start=rep(0, N_clean), t_stop=delays[['clean']], t_truncate=rep(Inf, N_clean),
  kappa=rep(1,3), n_parameters=3, 
  sigma_log=1/2, sigma_identity=1, sigma_logit=2/3
) %>% my_rdump(file='generalized-gamma-on-gamma-data.rdump')

generalized_gamma_on_gesgm_data <- list(
  n_obs=N_contaminated,
  t_start=rep(0, N_contaminated), t_stop=delays[['contaminated']], t_truncate=rep(Inf, N_contaminated),
  kappa=rep(1,3), n_parameters=3, 
  sigma_log=1/2, sigma_identity=1, sigma_logit=2/3
) %>% my_rdump(file='generalized-gamma-on-gesgm-data.rdump')

gesgm_on_gamma_data <- list(
  n_obs=N_clean,
  t_start=rep(0, N_clean), t_stop=delays[['clean']], t_truncate=rep(Inf, N_clean),
  kappa=c(10^-1, 10^-1, 1, .2, 1, 1, 1, 1), n_parameters=8, 
  sigma_log=1/2, sigma_identity=1, sigma_logit=2/3
) %>% my_rdump(file='gesgm-on-gamma-data.rdump')

gesgm_on_gesgm_data <- list(n_obs=N_contaminated,
  t_start=rep(0, N_contaminated), t_stop=delays[['contaminated']], t_truncate=rep(Inf, N_contaminated),
  kappa=c(1, 1, 1, .2, .1, .1, .1, 1), n_parameters=8, 
  sigma_log=1/2, sigma_identity=1, sigma_logit=2/3
) %>% my_rdump(file='gesgm-on-gesgm-data.rdump')


