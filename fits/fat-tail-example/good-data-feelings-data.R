
generate_delays <- function(N, alpha, beta, q, delta) {
  st_delays <- rgamma(n=N, shape=alpha, scale=beta)
  lt_delays <- ifelse(runif(n=N) < q, st_delays + rexp(n=N, rate=1/delta), st_delays)
  return(list(clean=st_delays, contaminated=lt_delays))
}


delay_data_plots <- function(delays) {
  require(dplyr); require(magrittr); require(ggplot2)

  data_q <- lapply(delays, quantile, probs=seq(from=.05, to=.95, by=0.01)) %>% data.frame
  pl_qq <- ggplot(data=data_q, aes(x=clean, y=contaminated)) + 
    geom_point() + geom_abline(slope=1, intercept=0)

  data <- data.frame(delays) %>% tidyr::gather(tail, delay)
  pl_histogram <- ggplot(data=data, aes(x=delay)) + 
    geom_histogram(bins=100) + 
    facet_wrap( ~ tail, ncol=1) + theme_minimal()
 
  return(list(qq=pl_qq, histogram=pl_histogram))
}

N <- 10^3; alpha <- 6; beta <- 1; q <- .05; delta=15

delays <- generate_delays(N, alpha, beta, q, delta)
delay_plots <- delay_data_plots(delays)
saveRDS(delays, file='delay-data.rds')
saveRDS(delay_plots, file='delay-data-plots.rds')



