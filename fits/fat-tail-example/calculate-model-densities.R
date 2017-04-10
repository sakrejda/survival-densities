# Function notes:
# 1) Input is parameter samples file.
# 2) Extracting parameters is model-specific
#   - This part could be skipped in a functionalized vesrion.
# 3) Calculating density with parameters is model-specific
#   - This calculation could pull parameters directly from the sample d.f.
# 4) Create a matrix n_samples x n_pts
# 5) Fill it row-wise looping over samples.
# 6) Apply on columns to get quantile, na.rm=TRUE estimates of the function (density
#    in this case).
# 7) Reduction to long format via tidy is standard.

library(magrittr); library(stannis)

pts <- seq(from=0, to=80, by=0.01)
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

gamma_on_gamma_densities <- matrix(data=NA, nrow=nrow(gamma_on_gamma_samples), ncol=length(pts),
  dimnames=list(iteration=gamma_on_gamma_samples[['iteration']], density=pts)
)
for ( i in 1:nrow(gamma_on_gamma_samples)) {
  alpha <- gamma_on_gamma_samples[['g_alpha']][i]
  beta <- gamma_on_gamma_samples[['g_beta']][i]
  gamma_on_gamma_densities[i,] <- dgamma(x=pts, shape=alpha, scale=beta)
}
gamma_on_gamma_densities_est <- apply(gamma_on_gamma_densities, 2, 
  quantile, na.rm=TRUE, probs=c(0.1, 0.5, 0.9)) %>% t %>% 
  data.frame(check.names=FALSE) %>% 
  mutate(density = as.numeric(rownames(.))) %>% 
  tidyr::gather(estimate, value, -density) 

# Gamma model on contaminated data:
gamma_on_gesgm_samples <- output[['gamma-on-gesgm']][['data']] %>% filter(stationary)

gamma_on_gesgm_densities <- matrix(data=NA, nrow=nrow(gamma_on_gesgm_samples), ncol=length(pts),
  dimnames=list(iteration=gamma_on_gesgm_samples[['iteration']], density=pts)
)
for ( i in 1:nrow(gamma_on_gesgm_samples)) {
  alpha <- gamma_on_gesgm_samples[['g_alpha']][i]
  beta <- gamma_on_gesgm_samples[['g_beta']][i]
  gamma_on_gesgm_densities[i,] <- dgamma(x=pts, shape=alpha, scale=beta)
}
gamma_on_gesgm_densities_est <- apply(gamma_on_gesgm_densities, 2, 
  quantile, na.rm=TRUE, probs=c(0.1, 0.5, 0.9)) %>% t %>% 
  data.frame(check.names=FALSE) %>% 
  mutate(density = as.numeric(rownames(.))) %>% 
  tidyr::gather(estimate, value, -density) 

# Generalized gamma model on clean data:
generalized_gamma_on_gamma_samples <- output[['generalized-gamma-on-gamma']][['data']] %>% 
  filter(stationary)

generalized_gamma_on_gamma_densities <- matrix(data=NA, nrow=nrow(generalized_gamma_on_gamma_samples), 
  ncol=length(pts), dimnames=list(iteration=generalized_gamma_on_gamma_samples[['iteration']], 
  density=pts)
)

for ( i in 1:nrow(generalized_gamma_on_gamma_samples)) {
  alpha <- generalized_gamma_on_gamma_samples[['g_alpha']][i]
  beta <- generalized_gamma_on_gamma_samples[['g_beta']][i]
  nu <- generalized_gamma_on_gamma_samples[['g_nu']][i]
  generalized_gamma_on_gamma_densities[i,] <- waitup:::generalized_gamma_pdf_1(x=pts, alpha=alpha, beta=beta, nu=nu)
}
generalized_gamma_on_gamma_densities_est <- apply(generalized_gamma_on_gamma_densities, 2, 
  quantile, na.rm=TRUE, probs=c(0.1, 0.5, 0.9)) %>% t %>% 
  data.frame(check.names=FALSE) %>% 
  mutate(density = as.numeric(rownames(.))) %>% 
  tidyr::gather(estimate, value, -density) 

# Generalized gamma model on contaminated data:
generalized_gamma_on_gesgm_samples <- output[['generalized-gamma-on-gesgm']][['data']] %>% 
  filter(stationary)

generalized_gamma_on_gesgm_densities <- matrix(data=NA, nrow=nrow(generalized_gamma_on_gesgm_samples), 
  ncol=length(pts), dimnames=list(iteration=generalized_gamma_on_gesgm_samples[['iteration']], 
  density=pts)
)
for ( i in 1:nrow(generalized_gamma_on_gesgm_samples)) {
  alpha <- generalized_gamma_on_gesgm_samples[['g_alpha']][i]
  beta <- generalized_gamma_on_gesgm_samples[['g_beta']][i]
  nu <- generalized_gamma_on_gesgm_samples[['g_nu']][i]
  generalized_gamma_on_gesgm_densities[i,] <- waitup:::generalized_gamma_pdf_1(x=pts, alpha=alpha, beta=beta, nu=nu) 
}
generalized_gamma_on_gesgm_densities_est <- apply(generalized_gamma_on_gesgm_densities, 2, 
  quantile, na.rm=TRUE, probs=c(0.1, 0.5, 0.9)) %>% t %>% 
  data.frame(check.names=FALSE) %>% 
  mutate(density = as.numeric(rownames(.))) %>% 
  tidyr::gather(estimate, value, -density) 

# Gamma exponential sum gamma mixture model on clean data:
gesgm_on_gamma_samples <- output[['gesgm-on-gamma']][['data']] %>% 
  filter(stationary)

gesgm_on_gamma_densities <- matrix(data=NA, nrow=nrow(gesgm_on_gamma_samples), 
  ncol=length(pts), dimnames=list(iteration=gesgm_on_gamma_samples[['iteration']], 
  density=pts)
)
for ( i in 1:nrow(gesgm_on_gamma_samples)) {
  q <- gesgm_on_gamma_samples[['g_q']][i]
  alpha <- gesgm_on_gamma_samples[['g_alpha']][i]
  beta <- gesgm_on_gamma_samples[['g_beta']][i]
  delta <- gesgm_on_gamma_samples[['g_delta']][i]
  gesgm_on_gamma_densities[i,] <- waitup:::gamma_exp_sum_gamma_mix_pdf_1(x=pts, q=q, alpha=alpha, beta=beta, delta=delta) 
}
gesgm_on_gamma_densities_est <- apply(gesgm_on_gamma_densities, 2, 
  quantile, na.rm=TRUE, probs=c(0.1, 0.5, 0.9)) %>% t %>% 
  data.frame(check.names=FALSE) %>% 
  mutate(density = as.numeric(rownames(.))) %>% 
  tidyr::gather(estimate, value, -density) 

# Gamma exponential sum gamma mixture model on contaminated data:
gesgm_on_gesgm_samples <- output[['gesgm-on-gesgm']][['data']] %>% 
  filter(stationary)

gesgm_on_gesgm_densities <- matrix(data=NA, nrow=nrow(gesgm_on_gesgm_samples), 
  ncol=length(pts), dimnames=list(iteration=gesgm_on_gesgm_samples[['iteration']], 
  density=pts)
)
for ( i in 1:nrow(gesgm_on_gesgm_samples)) {
  q <- gesgm_on_gesgm_samples[['g_q']][i]
  alpha <- gesgm_on_gesgm_samples[['g_alpha']][i]
  beta <- gesgm_on_gesgm_samples[['g_beta']][i]
  delta <- gesgm_on_gesgm_samples[['g_delta']][i]
  gesgm_on_gesgm_densities[i,] <- waitup:::gamma_exp_sum_gamma_mix_pdf_1(x=pts, q=q, alpha=alpha, beta=beta, delta=delta) 
}
gesgm_on_gesgm_densities_est <- apply(gesgm_on_gesgm_densities, 2, 
  quantile, na.rm=TRUE, probs=c(0.1, 0.5, 0.9)) %>% t %>% 
  data.frame(check.names=FALSE) %>% 
  mutate(density = as.numeric(rownames(.))) %>% 
  tidyr::gather(estimate, value, -density) 



densities <- list(
  gamma = list(
    clean = gamma_on_gamma_densities_est,
    contaminated = gamma_on_gesgm_densities_est
  ),
  generalized_gamma = list(
    clean = generalized_gamma_on_gamma_densities_est,
    contaminated = generalized_gamma_on_gesgm_densities_est
  ),
  gesgm = list(
    clean = gesgm_on_gamma_densities_est,
    contaminated = gesgm_on_gesgm_densities_est
  )
)




saveRDS(densities, file='model-densities.rds')





