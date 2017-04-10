library(rstan); library(dplyr); library(magrittr)
n_iterations <- 2000
n_warmup <- 1500
n_draws <- n_iterations - n_warmup

N <- 10^3
st_delays <- flexsurv::rgengamma.orig(n=N, shape=1, scale=1, k=3)
lt_delays <- ifelse(runif(n=N) < .05, st_delays + rexp(n=N, rate=1/5), st_delays)

delay_quantiles <- seq(from=.05, to=.95, by=0.01)
st_delays_q <- quantile(st_delays, probs=delay_quantiles)
lt_delays_q <- quantile(lt_delays, probs=delay_quantiles)

# Them quantiles are not on the 1:1 line.
plot(st_delays_q, lt_delays_q, xlim=c(0,5), ylim=c(0,5)) 
abline(a=0, b=1)

data <- data.frame(long=lt_delays, short=st_delays) %>%
  tidyr::gather(tail, delay)

pl_delays <- ggplot(data=data, aes(x=delay)) + 
  geom_histogram(bins=100) + 
  facet_wrap( ~ tail, ncol=1) + theme_minimal()
print(pl_delays)


compiled_gamma <- stan_model(file='../../models/stan-lang/full/gamma-p1.stan', 
  model_name='gamma-p1', save_dso=TRUE, obfuscate_model_name=FALSE)

compiled_gesgm <- stan_model(file='../../models/stan-lang/full/gamma-exp-sum-gamma-mix-p2.stan', 
  model_name='gamma-exp-sum-gamma-mix-p2', save_dso=TRUE, obfuscate_model_name=FALSE)

# First we use the regular gamma model on the gamma data set:
if (!exists('gamma_on_gamma')) {
  gamma_on_gamma <- sampling(compiled_gamma, data=list(n_obs=N, n_parameters=2,
    t_start=rep(0, N), t_stop=st_delays, t_truncate=rep(Inf, N),
    kappa=rep(1,2), sigma_log=1/2, sigma_identity=1, sigma_logit=2/3),
    chains=5, iter=n_iterations, warmup=n_warmup, cores=5,
    sample_file=paste0('gamma-on-gamma-out.'),
    diagnostic_file=paste0('gamma-on-gamma-diagnostic.')
  )
}

# Next we use the regular gamma model on the heavier-tailed dataset:
if (!exists('gamma_on_gesgm')) {
  gamma_on_gesgm <- sampling(compiled_gamma, data=list(n_obs=N, n_parameters=2,
    t_start=rep(0, N), t_stop=lt_delays, t_truncate=rep(Inf, N),
    kappa=rep(1,2), sigma_log=1/2, sigma_identity=1, sigma_logit=2/3),
    chains=5, iter=n_iterations, warmup=n_warmup, cores=5,
    sample_file=paste0('gamma-on-gesgm-out.'),
    diagnostic_file=paste0('gamma-on-gesgm-diagnostic.')
  )
}

# Lets do a visual comparison of the two fits:
g_g_alpha <- rstan::extract(gamma_on_gamma)[['alpha_0']] %>% exp
g_g_beta <- rstan::extract(gamma_on_gamma)[['beta_0']] %>% exp
pl_g_st_delays <- ggplot(data=data %>% filter(tail=='short'), 
  aes(x=delay)) + geom_histogram(aes(delay, ..density..), bins=100) +
  scale_x_continuous("delay", limits=c(0,8))
for (i in sample(x=(n_warmup+1):n_iterations, size=10)) 
  pl_g_st_delays <- pl_g_st_delays + stat_function(fun=dgamma, 
    alpha=.5, args=list(shape=g_g_alpha[i], scale=g_g_beta[i]))

g_gesgm_alpha <- rstan::extract(gamma_on_gesgm)[['alpha_0']] %>% exp
g_gesgm_beta <- rstan::extract(gamma_on_gesgm)[['beta_0']] %>% exp
pl_g_lt_delays <- ggplot(data=data %>% filter(tail=='long'), 
  aes(x=delay)) + geom_histogram(aes(delay, ..density..), bins=100) +
  scale_x_continuous("delay", limits=c(0,8))
for (i in sample(x=(n_warmup+1):n_iterations, size=10)) 
  pl_g_lt_delays <- pl_g_lt_delays + stat_function(fun=dgamma, 
    alpha=.5, args=list(shape=g_gesgm_alpha[i], scale=g_gesgm_beta[i]))

# How's it look:
pdf(file='gamma-model-st-lt-comparison.pdf', width=6, height=8)
gridExtra::grid.arrange(pl_g_st_delays, pl_g_lt_delays)
dev.off()

# That sugesgmests the gamma model is not doing so good with longer tailed
# data... but it don't look so bad... how bad is it? Let's look at a q-q plot:
g_st_delays_tq <- vector(mode='numeric', length=length(st_delays_q))
for (i in seq_along(delay_quantiles)) {
  g_st_delays_tq[i] <- qgamma(p=delay_quantiles[i],
    shape=mean(g_g_alpha), scale=mean(g_g_beta))
}

g_lt_delays_tq <- vector(mode='numeric', length=length(lt_delays_q))
for (i in seq_along(delay_quantiles)) {
  g_lt_delays_tq[i] <- qgamma(p=delay_quantiles[i],
    shape=mean(g_gesgm_alpha), scale=mean(g_gesgm_beta))
}

g_quantiles <- rbind(st_delays_q, g_st_delays_tq, lt_delays_q, g_lt_delays_tq) %>%
  data.frame(check.names=FALSE) %>% mutate(
    source=c(rep('short',2), rep('long',2)), 
    type=c('data','model','data','model')) %>%
  tidyr::gather(percentile, value, -source, -type) %>%
  tidyr::spread(type, value)

pl_g_st_lt_qq_comparison <- ggplot(data=g_quantiles, 
  aes(x=data, y=model, colour=percentile)
) + geom_point() +
    facet_wrap( ~ source, ncol=1) +
    geom_abline(slope=1, intercept=0) + 
    theme_minimal() + guides(colour=FALSE) +
    scale_x_continuous("delay (data quantiles)", limits=c(0,5)) +
    scale_y_continuous("delay (theoretical quantiles)", limits=c(0,5)) +
    coord_fixed()

# How's it look:
pdf(file='gamma-model-st-lt-qq-comparison.pdf', width=6, height=8)
print(pl_g_st_lt_qq_comparison)
dev.off()

# Dump these to screen also:
gridExtra::grid.arrange(pl_g_st_delays, pl_g_lt_delays)
print(pl_g_st_lt_qq_comparison)


# Now we use the generalized gamma model on the gamma data set:
gesgm_on_gamma_inits <- function() list(
  mu_0 = sqrt(3), sigma_0=sqrt(sqrt(3)), q_0=-5, delta_0=0)

if (!exists('gesgm_on_gamma')) {
  gesgm_on_gamma <- sampling(compiled_gesgm, data=list(n_obs=N,
    t_start=rep(0, N), t_stop=st_delays, t_truncate=rep(Inf, N),
    kappa=c(10^-1, 10^-1, 1, .2, 1, 1, 1, 1, 2), n_parameters=8, 
    sigma_log=1/2, sigma_identity=1, sigma_logit=2/3),
    chains=5, iter=n_iterations, warmup=n_warmup, cores=5,
    init=gesgm_on_gamma_inits,
    sample_file=paste0('gesgm-on-gamma-out.'),
    diagnostic_file=paste0('gesgm-on-gamma-diagnostic.')
  )
}


# Next we use the generalized gamma model on the heavier-tailed dataset:
gesgm_on_gesgm_inits <- function() list(
  mu_0 =log(sqrt(3)*sqrt(3)), sigma_0=log(sqrt(3)), q_0=-1, delta_0=log(5/1)-3)

if (!exists('gesgm_on_gesgm')) {
  gesgm_on_gesgm <- sampling(compiled_gesgm, data=list(n_obs=N,
    t_start=rep(0, N), t_stop=lt_delays, t_truncate=rep(Inf, N),
    kappa=c(1, 1, 1, .2, .001, .001, .01, .001, 4), n_parameters=8, 
    sigma_log=1/2, sigma_identity=1, sigma_logit=2/3),
    chains=5, iter=n_iterations, warmup=n_warmup, cores=5, 
    control=list(max_treedepth=14),
    init=gesgm_on_gesgm_inits,
    sample_file=paste0('gesgm-on-gesgm-out.'),
    diagnostic_file=paste0('gesgm-on-gesgm-diagnostic.')
  )
}

# Lets do a visual comparison of the two fits:
f <- function(x, q, alpha, beta, delta) sapply(x, waitup:::gamma_exp_sum_gamma_mix_pdf_1S, 
  q=q, alpha=alpha, beta=beta, delta=delta)

gesgm_g_alpha <- rstan::extract(gesgm_on_gamma)[['alpha_0']] %>% exp
gesgm_g_beta <- rstan::extract(gesgm_on_gamma)[['beta_0']] %>% exp
gesgm_g_delta <- rstan::extract(gesgm_on_gamma)[['delta_0']] %>% exp
gesgm_g_q <- rstan::extract(gesgm_on_gamma)[['q_0']] %>% plogis 
pl_gesgm_st_delays <- ggplot(data=data %>% filter(tail=='short'), 
  aes(x=delay)) + geom_histogram(aes(delay, ..density..), bins=100) +
  scale_x_continuous("delay", limits=c(0,8))
for (i in sample(x=(n_warmup+1):n_iterations, size=10)) 
  pl_gesgm_st_delays <- pl_gesgm_st_delays + stat_function(fun=f, alpha=.5, 
    args=list(q=gesgm_g_q[i], alpha=gesgm_g_delta[i], beta=gesgm_g_beta[i], delta=gesgm_g_beta[i]+gesgm_g_delta[i]))

gesgm_gesgm_alpha <- rstan::extract(gesgm_on_gesgm)[['alpha_0']] %>% exp
gesgm_gesgm_beta <- rstan::extract(gesgm_on_gesgm)[['beta_0']] %>% exp
gesgm_gesgm_delta <- rstan::extract(gesgm_on_gesgm)[['delta_0']] %>% exp
gesgm_gesgm_q <- rstan::extract(gesgm_on_gesgm)[['q_0']] %>% exp
pl_gesgm_lt_delays <- ggplot(data=data %>% filter(tail=='long'), 
  aes(x=delay)) + geom_histogram(aes(delay, ..density..), bins=100) +
  scale_x_continuous("delay", limits=c(0,8))
for (i in sample(x=(n_warmup+1):n_iterations, size=10)) 
  pl_gesgm_lt_delays <- pl_gesgm_lt_delays + stat_function(fun=flexsurv::dgengamma.orig, alpha=.5, 
    args=list(shape=gesgm_gesgm_delta[i], scale=gesgm_gesgm_beta[i], k=gesgm_gesgm_alpha[i]/gesgm_gesgm_delta[i]))

# How's it look:
pdf(file='gesgm-model-st-lt-comparison.pdf', width=6, height=8)
gridExtra::grid.arrange(pl_gesgm_st_delays, pl_gesgm_lt_delays)
dev.off()

# That sugesgmests the gamma model is not doing so good with longer tailed
# data... but it don't look so bad... how bad is it? Let's look at a q-q plot:
gesgm_st_delays_tq <- vector(mode='numeric', length=length(st_delays_q))
for (i in seq_along(delay_quantiles)) {
  gesgm_st_delays_tq[i] <- flexsurv::qgengamma.orig(p=delay_quantiles[i],
    shape=mean(gesgm_g_delta), scale=mean(gesgm_g_beta), k=mean(gesgm_g_alpha/gesgm_g_delta))
}

gesgm_lt_delays_tq <- vector(mode='numeric', length=length(lt_delays_q))
for (i in seq_along(delay_quantiles)) {
  gesgm_lt_delays_tq[i] <- flexsurv::qgengamma.orig(p=delay_quantiles[i],
    shape=mean(gesgm_gesgm_delta), scale=mean(gesgm_gesgm_beta), k=mean(gesgm_gesgm_alpha/gesgm_gesgm_delta))
}

gesgm_quantiles <- rbind(st_delays_q, gesgm_st_delays_tq, lt_delays_q, gesgm_lt_delays_tq) %>%
  data.frame(check.names=FALSE) %>% mutate(
    source=c(rep('short',2), rep('long',2)), 
    type=c('data','model','data','model')) %>%
  tidyr::gather(percentile, value, -source, -type) %>%
  tidyr::spread(type, value)

pl_gesgm_st_lt_qq_comparison <- ggplot(data=gesgm_quantiles, 
  aes(x=data, y=model, colour=percentile)
) + geom_point() +
    facet_wrap( ~ source, ncol=1) +
    geom_abline(slope=1, intercept=0) + 
    theme_minimal() + guides(colour=FALSE) +
    scale_x_continuous("delay (data quantiles)", limits=c(0,2.5)) +
    scale_y_continuous("delay (theoretical quantiles)", limits=c(0,2.5)) +
    coord_fixed()

# How's it look:
pdf(file='gesgm-model-st-lt-qq-comparison.pdf', width=6, height=8)
print(pl_gesgm_st_lt_qq_comparison)
dev.off()

# Dump these to screen also:
gridExtra::grid.arrange(pl_gesgm_st_delays, pl_gesgm_lt_delays)
print(pl_gesgm_st_lt_qq_comparison)


# Four-way comarison by densities
pdf(file='model-st-lt-density-comparison.pdf', width=6, height=8)
gridExtra::grid.arrange(pl_g_st_delays, pl_g_lt_delays,
  pl_gesgm_st_delays, pl_gesgm_lt_delays, ncol=2)
dev.off()


all_quantiles <- rbind(
  g_quantiles %>% mutate(fit='gamma'),
  gesgm_quantiles %>% mutate(fit='gesgm')
)

pl_st_lt_qq_comparison <- ggplot(
  data=all_quantiles,
  aes(x=data, y=model, colour=fit, group=paste(fit, percentile))
) + geom_point() +
    facet_wrap( ~ source, ncol=2) +
    geom_abline(slope=1, intercept=0) + 
    theme_minimal() + guides(colour=FALSE) +
    scale_x_continuous("delay (data quantiles)", limits=c(0,2.5)) +
    scale_y_continuous("delay (theoretical quantiles)", limits=c(0,2.5)) +
    coord_fixed()

# Four-way comarison by q-q plots
pdf(file='model-st-lt-qq-comparison.pdf', width=6, height=8)
print(pl_st_lt_qq_comparison)
dev.off()





# Exploratory plots for gesgm-on-gamma:
start_iter <- 0
library(ggplot2)
output <- dir(pattern='gesgm-on-gamma-out')
o <- lapply(output, function(x) {
  a <- read.csv(file=x, stringsAsFactors=FALSE, comment.char='#')
  a$file=x
  a$iteration=1:nrow(a)
  return(a)}
)
om <- do.call(what=rbind, args=o)
ggplot(data=om[om$iteration>start_iter,], aes(x=iteration, y=log10(stepsize__), alpha=iteration)) + geom_point(size=0.001) + facet_wrap( ~ file)
ggplot(data=om[om$iteration>start_iter,], aes(x=iteration, y=mu_0, alpha=iteration)) + geom_point(size=0.001) + facet_wrap( ~ file)
ggplot(data=om[om$iteration>start_iter,], aes(x=iteration, y=sigma_0, alpha=iteration)) + geom_point(size=0.001) + facet_wrap( ~ file)
ggplot(data=om[om$iteration>start_iter,], aes(x=iteration, y=delta_0, alpha=iteration)) + geom_point(size=0.001) + facet_wrap( ~ file)
ggplot(data=om[om$iteration>start_iter,], aes(x=iteration, y=q_0, alpha=iteration)) + geom_point(size=0.001) + facet_wrap( ~ file)


ggplot(data=om, aes(x=alpha_0, y=beta_0)) + geom_point() + facet_wrap( ~ file)
ggplot(data=om, aes(x=iteration, y=alpha_0)) + geom_point() + facet_wrap( ~ file)
ggplot(data=om[om$iteration>start_iter,], aes(x=iteration, y=alpha_0)) + geom_point() + facet_wrap( ~ file)
ggplot(data=om[om$iteration>start_iter,], aes(x=alpha_0, y=beta_0)) + geom_point() + facet_wrap( ~ file)
ggplot(data=om[om$iteration>start_iter,], aes(x=alpha_0, y=sigma_0)) + geom_point() + facet_wrap( ~ file)

ggplot(data=om, aes(x=alpha_0, y=beta_0)) + geom_point() + facet_wrap( ~ file)
ggplot(data=om[om$iteration>start_iter,], aes(x=iteration, y=(stepsize__), alpha=iteration)) + geom_point(size=0.001) + facet_wrap( ~ file)
ggplot(data=om[om$iteration>start_iter,], aes(x=iteration, y=(accept_stat__), alpha=iteration)) + geom_point(size=0.001) + facet_wrap( ~ file)
ggplot(data=om[om$iteration>start_iter,], aes(x=iteration, y=(accept_stat__), alpha=iteration, colour=divergent__)) + geom_point(size=0.001) + facet_wrap( ~ file)
ggplot(data=om[om$iteration>start_iter,], aes(x=iteration, y=(accept_stat__), colour=divergent__)) + geom_point(size=0.001) + facet_wrap( ~ file)
ggplot(data=om[om$iteration>start_iter,], aes(x=iteration, y=(accept_stat__), colour=(1==divergent__))) + geom_point(size=0.001) + facet_wrap( ~ file)
ggplot(data=om[om$iteration>start_iter,], aes(x=iteration, y=alpha_0, colour=accept_stat__)) + geom_point(size=0.001) + facet_wrap( ~ file)
ggplot(data=om[om$iteration>start_iter,], aes(x=iteration, y=q_0, colour=accept_stat__)) + geom_point(size=0.001) + facet_wrap( ~ file)
ggplot(data=om[om$iteration>start_iter,], aes(x=iteration, y=beta_0, colour=accept_stat__)) + geom_point(size=0.001) + facet_wrap( ~ file)
ggplot(data=om[om$iteration>start_iter,], aes(x=alpha_0, y=beta_0, colour=accept_stat__)) + geom_point(size=0.001) + facet_wrap( ~ file)
ggplot(data=om[om$iteration>start_iter,], aes(x=mu_0, y=sigma_0, colour=accept_stat__)) + geom_point(size=0.001) + facet_wrap( ~ file)
ggplot(data=om[om$iteration>start_iter,], aes(x=mu_0, y=delta_0, colour=accept_stat__)) + geom_point(size=0.001) + facet_wrap( ~ file)
ggplot(data=om[om$iteration>start_iter,], aes(x=sigma_0, y=delta_0, colour=accept_stat__)) + geom_point(size=0.001) + facet_wrap( ~ file)



# Exploratory plots for gesgm-on-gesgm:
start_iter <- 0
library(ggplot2)
output <- dir(pattern='gesgm-on-gesgm-out')
o <- lapply(output, function(x) {
  a <- read.csv(file=x, stringsAsFactors=FALSE, comment.char='#')
  a$file=x
  a$iteration=1:nrow(a)
  return(a)}
)
om <- do.call(what=rbind, args=o)
ggplot(data=om[om$iteration>start_iter,], aes(x=iteration, y=log10(stepsize__), alpha=iteration)) + geom_point(size=0.001) + facet_wrap( ~ file)
ggplot(data=om[om$iteration>start_iter,], aes(x=iteration, y=mu_0, alpha=iteration)) + geom_point(size=0.001) + facet_wrap( ~ file)
ggplot(data=om[om$iteration>start_iter,], aes(x=iteration, y=sigma_0, alpha=iteration)) + geom_point(size=0.001) + facet_wrap( ~ file)
ggplot(data=om[om$iteration>start_iter,], aes(x=iteration, y=delta_0, alpha=iteration)) + geom_point(size=0.001) + facet_wrap( ~ file)
ggplot(data=om[om$iteration>start_iter,], aes(x=iteration, y=q_0, alpha=iteration)) + geom_point(size=0.001) + facet_wrap( ~ file)



ggplot(data=om, aes(x=alpha_0, y=beta_0)) + geom_point() + facet_wrap( ~ file)
ggplot(data=om, aes(x=iteration, y=alpha_0)) + geom_point() + facet_wrap( ~ file)
ggplot(data=om[om$iteration>start_iter,], aes(x=iteration, y=alpha_0)) + geom_point() + facet_wrap( ~ file)
ggplot(data=om[om$iteration>start_iter,], aes(x=alpha_0, y=beta_0)) + geom_point() + facet_wrap( ~ file)
ggplot(data=om[om$iteration>start_iter,], aes(x=alpha_0, y=sigma_0)) + geom_point() + facet_wrap( ~ file)

ggplot(data=om, aes(x=alpha_0, y=beta_0)) + geom_point() + facet_wrap( ~ file)
ggplot(data=om[om$iteration>start_iter,], aes(x=iteration, y=(stepsize__), alpha=iteration)) + geom_point(size=0.001) + facet_wrap( ~ file)
ggplot(data=om[om$iteration>start_iter,], aes(x=iteration, y=(accept_stat__), alpha=iteration)) + geom_point(size=0.001) + facet_wrap( ~ file)
ggplot(data=om[om$iteration>start_iter,], aes(x=iteration, y=(accept_stat__), alpha=iteration, colour=divergent__)) + geom_point(size=0.001) + facet_wrap( ~ file)
ggplot(data=om[om$iteration>start_iter,], aes(x=iteration, y=(accept_stat__), colour=divergent__)) + geom_point(size=0.001) + facet_wrap( ~ file)
ggplot(data=om[om$iteration>start_iter,], aes(x=iteration, y=(accept_stat__), colour=(1==divergent__))) + geom_point(size=0.001) + facet_wrap( ~ file)
ggplot(data=om[om$iteration>start_iter,], aes(x=iteration, y=alpha_0, colour=accept_stat__)) + geom_point(size=0.001) + facet_wrap( ~ file)
ggplot(data=om[om$iteration>start_iter,], aes(x=iteration, y=q_0, colour=accept_stat__)) + geom_point(size=0.001) + facet_wrap( ~ file)
ggplot(data=om[om$iteration>start_iter,], aes(x=iteration, y=beta_0, colour=accept_stat__)) + geom_point(size=0.001) + facet_wrap( ~ file)
ggplot(data=om[om$iteration>start_iter,], aes(x=alpha_0, y=beta_0, colour=accept_stat__)) + geom_point(size=0.001) + facet_wrap( ~ file)
ggplot(data=om[om$iteration>start_iter,], aes(x=mu_0, y=sigma_0, colour=accept_stat__)) + geom_point(size=0.001) + facet_wrap( ~ file)
ggplot(data=om[om$iteration>start_iter,], aes(x=mu_0, y=delta_0, colour=accept_stat__)) + geom_point(size=0.001) + facet_wrap( ~ file)
ggplot(data=om[om$iteration>start_iter,], aes(x=sigma_0, y=delta_0, colour=accept_stat__)) + geom_point(size=0.001) + facet_wrap( ~ file)

