library(magrittr); set.seed(42)

generate_delays <- function(N, alpha, beta, q, delta) {
  st_delays <- rgamma(n=N, shape=alpha, scale=beta)
  lt_delays <- ifelse(runif(n=N) < q, st_delays + rexp(n=N, rate=1/delta), st_delays)
  return(list(clean=st_delays, contaminated=lt_delays))
}


N <- 10^3; 

theta <- list(
  alpha = c(1, 1.5, 2, 10, 20),
  beta = c(0.1, 1, 1.5, 10, 20),
  q = c(0.05, 0.1, 0.5),
  delta = c(5, 10, 50, 100)
)

theta_grid <- do.call(what=expand.grid, args=theta)

delays <- list()
for (i in 1:nrow(theta_grid)) {
  delays[[i]] <- generate_delays(N=N, alpha=theta_grid[['alpha']],
    beta=theta_grid[['beta']], q=theta_grid[['q']], delta=theta_grid[['delta']])
}

saveRDS(delays, file='delay-data.rds')
saveRDS(theta_grid, file='theta.rds')



