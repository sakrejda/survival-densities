library(magrittr)

qps <- seq(from=0.05, to=0.95, by=0.01)

N <- 10^3; alpha <- 3; beta <- 1; q <- .05; delta=5

delay_data <- readRDS('delay-data.rds')
data_q <- lapply(delay_data, quantile, probs=qps) %>% data.frame

clean_theoretical_q <- qgamma(p=qps, shape=alpha, scale=beta)
contaminated_theoretical_q <- sapply(qps, function(p) {
  f <-function(x) waitup:::gamma_exp_sum_gamma_mix_cdf_1(
    x=x, q=q, alpha=alpha, beta=beta, delta=delta) - p
  q <- uniroot(f, interval=c(0,100), f.lower=-p, f.upper=1-p,
    check.conv=TRUE, maxiter=10^4)
  return(q$root)
})


quantiles <- list(
  data = list(
    clean = data.frame(quantile=qps, value=data_q[['clean']], data='clean', stringsAsFactors=FALSE),
    contaminated = data.frame(quantile=qps, value=data_q[['contaminated']], data='contaminated', stringsAsFactors=FALSE)
  ), 
  theoretical = list(
    clean = data.frame(quantile=qps, value=clean_theoretical_q, data='clean', stringsAsFactors=FALSE),
    contaminated = data.frame(quantile=qps, value=contaminated_theoretical_q, data='contaminated', stringsAsFactors=FALSE)
  )
) %>% lapply(function(x) do.call(what=rbind, args=x))

saveRDS(quantiles, file='data-quantiles.rds')





