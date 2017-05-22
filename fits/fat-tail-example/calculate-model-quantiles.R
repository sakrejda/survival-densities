library(magrittr); library(stannis)

qps <- seq(from=0.05, to=0.95, by=0.01)
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

gamma_on_gamma_quantiles <- matrix(data=NA, nrow=nrow(gamma_on_gamma_samples), ncol=length(qps),
  dimnames=list(iteration=gamma_on_gamma_samples[['iteration']], quantile=qps)
)
for ( i in 1:nrow(gamma_on_gamma_samples)) {
  alpha <- gamma_on_gamma_samples[['g_alpha']][i]
  beta <- gamma_on_gamma_samples[['g_beta']][i]
  gamma_on_gamma_quantiles[i,] <- qgamma(p=qps, shape=alpha, scale=beta)
}
gamma_on_gamma_quantiles_est <- apply(gamma_on_gamma_quantiles, 2, 
  quantile, probs=c(0.1, 0.5, 0.9)) %>% t %>% 
  data.frame(check.names=FALSE) %>% 
  mutate(quantile = as.numeric(rownames(.))) %>% 
  tidyr::gather(estimate, value, -quantile) 

# Gamma model on contaminated data:
gamma_on_gesgm_samples <- output[['gamma-on-gesgm']][['data']] %>% filter(stationary)

gamma_on_gesgm_quantiles <- matrix(data=NA, nrow=nrow(gamma_on_gesgm_samples), ncol=length(qps),
  dimnames=list(iteration=gamma_on_gesgm_samples[['iteration']], quantile=qps)
)
for ( i in 1:nrow(gamma_on_gesgm_samples)) {
  alpha <- gamma_on_gesgm_samples[['g_alpha']][i]
  beta <- gamma_on_gesgm_samples[['g_beta']][i]
  gamma_on_gesgm_quantiles[i,] <- qgamma(p=qps, shape=alpha, scale=beta)
}
gamma_on_gesgm_quantiles_est <- apply(gamma_on_gesgm_quantiles, 2, 
  quantile, probs=c(0.1, 0.5, 0.9)) %>% t %>% 
  data.frame(check.names=FALSE) %>% 
  mutate(quantile = as.numeric(rownames(.))) %>% 
  tidyr::gather(estimate, value, -quantile) 

# Generalized gamma model on clean data:
generalized_gamma_on_gamma_samples <- output[['generalized-gamma-on-gamma']][['data']] %>% 
  filter(stationary)

generalized_gamma_on_gamma_quantiles <- matrix(data=NA, nrow=nrow(generalized_gamma_on_gamma_samples), 
  ncol=length(qps), dimnames=list(iteration=generalized_gamma_on_gamma_samples[['iteration']], 
  quantile=qps)
)

generic_qdf <- function(cdf, domain=c(-Inf, Inf)) {
  cdf <- cdf
  domain <- domain
  f <- function(p, ...) {
    g <- function(q, ...) cdf(q, ...) - p
    q <- uniroot(g, interval=domain, f.lower=-p, f.upper=1-p,
                 check.conv=TRUE, maxiter=10^4, ...)
    return(q$root)
  }
  return(f)
}

for ( i in 1:nrow(generalized_gamma_on_gamma_samples)) {
  alpha <- generalized_gamma_on_gamma_samples[['g_alpha']][i]
  beta <- generalized_gamma_on_gamma_samples[['g_beta']][i]
  nu <- generalized_gamma_on_gamma_samples[['g_nu']][i]
  generalized_gamma_on_gamma_quantiles[i,] <- sapply(qps, function(p, alpha, beta, nu) {
    f <- function(x) waitup:::generalized_gamma_cdf_1(x=x, alpha=alpha, beta=beta, nu=nu) - p
    q <- uniroot(f, interval=c(0,50), f.lower=-p, f.upper=1-p,
      check.conv=TRUE, maxiter=10^4)
    return(q$root)
  }, alpha=alpha, beta=beta, nu=nu)
}
generalized_gamma_on_gamma_quantiles_est <- apply(generalized_gamma_on_gamma_quantiles, 2, 
  quantile, probs=c(0.1, 0.5, 0.9)) %>% t %>% 
  data.frame(check.names=FALSE) %>% 
  mutate(quantile = as.numeric(rownames(.))) %>% 
  tidyr::gather(estimate, value, -quantile) 

# Generalized gamma model on contaminated data:
generalized_gamma_on_gesgm_samples <- output[['generalized-gamma-on-gesgm']][['data']] %>% 
  filter(stationary)

generalized_gamma_on_gesgm_quantiles <- matrix(data=NA, nrow=nrow(generalized_gamma_on_gesgm_samples), 
  ncol=length(qps), dimnames=list(iteration=generalized_gamma_on_gesgm_samples[['iteration']], 
  quantile=qps)
)
for ( i in 1:nrow(generalized_gamma_on_gesgm_samples)) {
  alpha <- generalized_gamma_on_gesgm_samples[['g_alpha']][i]
  beta <- generalized_gamma_on_gesgm_samples[['g_beta']][i]
  nu <- generalized_gamma_on_gesgm_samples[['g_nu']][i]
  generalized_gamma_on_gesgm_quantiles[i,] <- sapply(qps, function(p, alpha, beta, nu) {
    f <- function(x) waitup:::generalized_gamma_cdf_1(x=x, alpha=alpha, beta=beta, nu=nu) - p
    q <- uniroot(f, interval=c(0,50), f.lower=-p, f.upper=1-p,
      check.conv=TRUE, maxiter=10^4)
    return(q$root)
  }, alpha=alpha, beta=beta, nu=nu)
}
generalized_gamma_on_gesgm_quantiles_est <- apply(generalized_gamma_on_gesgm_quantiles, 2, 
  quantile, probs=c(0.1, 0.5, 0.9)) %>% t %>% 
  data.frame(check.names=FALSE) %>% 
  mutate(quantile = as.numeric(rownames(.))) %>% 
  tidyr::gather(estimate, value, -quantile) 

# Gamma exponential sum gamma mixture model on clean data:
gesgm_on_gamma_samples <- output[['gesgm-on-gamma']][['data']] %>% 
  filter(stationary)

gesgm_on_gamma_quantiles <- matrix(data=NA, nrow=nrow(gesgm_on_gamma_samples), 
  ncol=length(qps), dimnames=list(iteration=gesgm_on_gamma_samples[['iteration']], 
  quantile=qps)
)
for ( i in 1:nrow(gesgm_on_gamma_samples)) {
  q <- gesgm_on_gamma_samples[['g_q']][i]
  alpha <- gesgm_on_gamma_samples[['g_alpha']][i]
  beta <- gesgm_on_gamma_samples[['g_beta']][i]
  delta <- gesgm_on_gamma_samples[['g_delta']][i]
  gesgm_on_gamma_quantiles[i,] <- sapply(qps, function(p, q, alpha, beta, delta) {
    f <- function(x) waitup:::gamma_exp_sum_gamma_mix_cdf_1(x=x, q=q, alpha=alpha, beta=beta, delta=delta) - p
    q <- uniroot(f, interval=c(0,50), f.lower=-p, f.upper=1-p,
      check.conv=TRUE, maxiter=10^4)
    return(q$root)
  }, q=q, alpha=alpha, beta=beta, delta=delta)
}
gesgm_on_gamma_quantiles_est <- apply(gesgm_on_gamma_quantiles, 2, 
  quantile, probs=c(0.1, 0.5, 0.9)) %>% t %>% 
  data.frame(check.names=FALSE) %>% 
  mutate(quantile = as.numeric(rownames(.))) %>% 
  tidyr::gather(estimate, value, -quantile) 

# Gamma exponential sum gamma mixture model on contaminated data:
gesgm_on_gesgm_samples <- output[['gesgm-on-gesgm']][['data']] %>% 
  filter(stationary)

gesgm_on_gesgm_quantiles <- matrix(data=NA, nrow=nrow(gesgm_on_gesgm_samples), 
  ncol=length(qps), dimnames=list(iteration=gesgm_on_gesgm_samples[['iteration']], 
  quantile=qps)
)
for ( i in 1:nrow(gesgm_on_gesgm_samples)) {
  q <- gesgm_on_gesgm_samples[['g_q']][i]
  alpha <- gesgm_on_gesgm_samples[['g_alpha']][i]
  beta <- gesgm_on_gesgm_samples[['g_beta']][i]
  delta <- gesgm_on_gesgm_samples[['g_delta']][i]
  gesgm_on_gesgm_quantiles[i,] <- sapply(qps, function(p, q, alpha, beta, delta) {
    f <- function(x) waitup:::gamma_exp_sum_gamma_mix_cdf_1(x=x, q=q, alpha=alpha, beta=beta, delta=delta) - p
    q <- uniroot(f, interval=c(0,50), f.lower=-p, f.upper=1-p,
      check.conv=TRUE, maxiter=10^4)
    return(q$root)
  }, q=q, alpha=alpha, beta=beta, delta=delta)
}
gesgm_on_gesgm_quantiles_est <- apply(gesgm_on_gesgm_quantiles, 2, 
  quantile, probs=c(0.1, 0.5, 0.9)) %>% t %>% 
  data.frame(check.names=FALSE) %>% 
  mutate(quantile = as.numeric(rownames(.))) %>% 
  tidyr::gather(estimate, value, -quantile) 



quantiles <- list(
  gamma = list(
    clean = gamma_on_gamma_quantiles_est,
    contaminated = gamma_on_gesgm_quantiles_est
  ),
  generalized_gamma = list(
    clean = generalized_gamma_on_gamma_quantiles_est,
    contaminated = generalized_gamma_on_gesgm_quantiles_est
  ),
  gesgm = list(
    clean = gesgm_on_gamma_quantiles_est,
    contaminated = gesgm_on_gesgm_quantiles_est
  )
)




saveRDS(quantiles, file='model-quantiles.rds')





