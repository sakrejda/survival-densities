library(magrittr); set.seed(42)

generate_delays <- function(N, alpha, beta, q, delta) {
  st_delays <- rgamma(n=N, shape=alpha, scale=beta)
  lt_delays <- ifelse(runif(n=N) < q, st_delays + rexp(n=N, rate=1/delta), st_delays)
  return(list(gamma=st_delays, gesgm=lt_delays))
}


N <- 10^3; 

theta <- list(
  alpha = c(1, 1.5, 8, 20),
  beta = 1,
  q = c(0.05, 0.1, 0.5),
  delta = 15,
  replicate = 1:100
)

theta_grid <- do.call(what=expand.grid, args=theta)

delays <- list()
for (i in 1:nrow(theta_grid)) {
  delays[[i]] <- generate_delays(N=N, alpha=theta_grid[i,'alpha'],
    beta=theta_grid[i,'beta'], q=theta_grid[i,'q'], delta=theta_grid[i,'delta'])
}

saveRDS(delays, file='delays.rds')
saveRDS(theta_grid, file='theta.rds')



