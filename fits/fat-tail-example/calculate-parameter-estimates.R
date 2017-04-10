library(magrittr); library(stannis)

qps <- seq(from=0.05, to=0.95, by=0.01)
est_probs <- c(0.1, 0.5, 0.9)
N <- 10^3; alpha <- 3; beta <- 1; q <- .05; delta=5

run_data <- readRDS(file='run-data.rds')
output_patterns <- '^' %>% 
  paste0(gsub(pattern='X', replacement='[0-9][0-9]', x=run_data$output_files)) %>%
  paste0('$')

diagnostic_patterns <-  gsub(pattern='X', replacement='[0-9][0-9]', x=run_data$diagnostic_files)

output <- lapply(output_patterns, function(x) read_stan_set(root='.', pattern=x)) %>%
  `names<-`(run_data[['run_names']])

# Gamma model on clean data:
gamma_on_gamma_samples <- output[['gamma-on-gamma']][['data']] %>% filter(stationary)

gamma_on_gamma_alpha_est <- gamma_on_gamma_samples[['g_alpha']] %>% quantile(probs=est_probs)
gamma_on_gamma_beta_est <- gamma_on_gamma_samples[['g_beta']] %>% quantile(probs=est_probs)

# Gamma model on contaminated data:
gamma_on_gesgm_samples <- output[['gamma-on-gesgm']][['data']] %>% filter(stationary)

gamma_on_gesgm_alpha_est <- gamma_on_gesgm_samples[['g_alpha']] %>% quantile(probs=est_probs)
gamma_on_gesgm_beta_est <- gamma_on_gesgm_samples[['g_beta']] %>% quantile(probs=est_probs)


# Generalized gamma model on clean data:
generalized_gamma_on_gamma_samples <- output[['generalized-gamma-on-gamma']][['data']] %>% 
  filter(stationary)

generalized_gamma_on_gamma_k_est <- generalized_gamma_on_gamma_samples[['g_k']] %>% quantile(probs=est_probs)
generalized_gamma_on_gamma_mu_est <- generalized_gamma_on_gamma_samples[['g_mu']] %>% quantile(probs=est_probs)
generalized_gamma_on_gamma_sigma_est <- generalized_gamma_on_gamma_samples[['g_sigma']] %>% quantile(probs=est_probs)

# Generalized gamma model on contaminated data:
generalized_gamma_on_gesgm_samples <- output[['generalized-gamma-on-gesgm']][['data']] %>% 
  filter(stationary)

generalized_gamma_on_gesgm_k_est <- generalized_gamma_on_gesgm_samples[['g_k']] %>% quantile(probs=est_probs)
generalized_gamma_on_gesgm_mu_est <- generalized_gamma_on_gesgm_samples[['g_mu']] %>% quantile(probs=est_probs)
generalized_gamma_on_gesgm_sigma_est <- generalized_gamma_on_gesgm_samples[['g_sigma']] %>% quantile(probs=est_probs)

# Gamma exponential sum gamma mixture model on clean data:
gesgm_on_gamma_samples <- output[['gesgm-on-gamma']][['data']] %>% 
  filter(stationary)

gesgm_on_gamma_alpha_est <- gesgm_on_gamma_samples[['g_alpha']] %>% quantile(probs=est_probs)
gesgm_on_gamma_beta_est <- gesgm_on_gamma_samples[['g_beta']] %>% quantile(probs=est_probs)
gesgm_on_gamma_delta_est <- gesgm_on_gamma_samples[['g_delta']] %>% quantile(probs=est_probs)

# Gamma exponential sum gamma mixture model on contaminated data:
gesgm_on_gesgm_samples <- output[['gesgm-on-gesgm']][['data']] %>% 
  filter(stationary)

gesgm_on_gesgm_alpha_est <- gesgm_on_gesgm_samples[['g_alpha']] %>% quantile(probs=est_probs)
gesgm_on_gesgm_beta_est <- gesgm_on_gesgm_samples[['g_beta']] %>% quantile(probs=est_probs)
gesgm_on_gesgm_delta_est <- gesgm_on_gesgm_samples[['g_delta']] %>% quantile(probs=est_probs)


parameter_estimates <- list(
  gamma_data = list(
    gamma_model = list(
      alpha=gamma_on_gamma_alpha_est, 
      beta=gamma_on_gamma_beta_est),
    gesgm_model = list(
      alpha=gesgm_on_gamma_alpha_est, 
      beta=gesgm_on_gamma_beta_est, 
      delta=gesgm_on_gamma_delta_est),
    generalized_gamma_model = list(
      k=generalized_gamma_on_gamma_k_est, 
      mu=generalized_gamma_on_gamma_mu_est, 
      sigma=generalized_gamma_on_gamma_sigma_est)
  ), 
  gesgm_data = list(
    gamma_model = list(
      alpha=gamma_on_gesgm_alpha_est, 
      beta=gamma_on_gesgm_beta_est),
    gesgm_model = list(
      alpha=gesgm_on_gesgm_alpha_est, 
      beta=gesgm_on_gesgm_beta_est, 
      delta=gesgm_on_gesgm_delta_est),
    generalized_gamma_model = list(
      k=generalized_gamma_on_gesgm_k_est, 
      mu=generalized_gamma_on_gesgm_mu_est, 
      sigma=generalized_gamma_on_gesgm_sigma_est)
  )
)

saveRDS(parameter_estimates, file='parameter-estimates.rds')




 



