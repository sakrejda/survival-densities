library(rstan); library(dplyr); library(magrittr)

N <- 10^3
st_delays <- flexsurv::rgengamma.orig(n=N, shape=3, scale=1, k=3/3)
st_delays <- st_delays/median(st_delays)
#lt_delays <- ifelse(runif(n=N) < .05, st_delays + runif(n=N, min=0, max=20), st_delays)
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


compiled_gamma <- stan_model(file='../../models/stan-lang/full-fix/gamma-p1.stan', 
  model_name='gamma-p1', save_dso=TRUE, obfuscate_model_name=FALSE)

compiled_generalized_gamma <- stan_model(file='../../models/stan-lang/full-fix/generalized-gamma-p1.stan', 
  model_name='generalized-gamma-p1', save_dso=TRUE, obfuscate_model_name=FALSE)

compiled_gamma_exp_sum_gamma_mix <- stan_model(file='../../models/stan-lang/full-fix/gamma-exp-sum-gamma-mix-p1.stan', 
  model_name='gamma-exp-sum-gamma-mix-p1', save_dso=TRUE, obfuscate_model_name=FALSE)

# First we use the regular gamma model on the gamma data set:
if (!exists('gamma_on_gamma')) {
  gamma_on_gamma <- sampling(compiled_gamma, data=list(n_obs=N, n_parameters=2,
    t_start=rep(0, N), t_stop=st_delays, t_truncate=rep(Inf, N),
    kappa=rep(1,2), sigma_log=1/2, sigma_identity=1, sigma_logit=2/3),
    chains=5, iter=800, sample_file=paste0('gamma-on-gamma-out.'),
    diagnostic_file=paste0('gamma-on-gamma-diagnostic.')
  )
}

# Next we use the regular gamma model on the heavier-tailed dataset:
if (!exists('gamma_on_generalized_gamma')) {
  gamma_on_generalized_gamma <- sampling(compiled_gamma, data=list(n_obs=N, n_parameters=2,
    t_start=rep(0, N), t_stop=lt_delays, t_truncate=rep(Inf, N),
    kappa=rep(1,2), sigma_log=1/2, sigma_identity=1, sigma_logit=2/3),
    chains=5, iter=800, sample_file=paste0('gamma-on-generalized-gamma-out.'),
    diagnostic_file=paste0('gamma-on-generalized-gamma-diagnostic.')
  )
}

# Lets do a visual comparison of the two fits:
g_g_alpha <- rstan::extract(gamma_on_gamma)[['alpha_0']] %>% exp
g_g_beta <- rstan::extract(gamma_on_gamma)[['beta_0']] %>% exp
pl_g_st_delays <- ggplot(data=data %>% filter(tail=='short'), 
  aes(x=delay)) + geom_histogram(aes(delay, ..density..), bins=100) +
  scale_x_continuous("delay", limits=c(0,8))
for (i in sample(x=1001:2000, size=10)) 
  pl_g_st_delays <- pl_g_st_delays + stat_function(fun=dgamma, 
    alpha=.5, args=list(shape=g_g_alpha[i], scale=g_g_beta[i]))

g_gg_alpha <- rstan::extract(gamma_on_generalized_gamma)[['alpha_0']] %>% exp
g_gg_beta <- rstan::extract(gamma_on_generalized_gamma)[['beta_0']] %>% exp
pl_g_lt_delays <- ggplot(data=data %>% filter(tail=='long'), 
  aes(x=delay)) + geom_histogram(aes(delay, ..density..), bins=100) +
  scale_x_continuous("delay", limits=c(0,8))
for (i in sample(x=1001:2000, size=10)) 
  pl_g_lt_delays <- pl_g_lt_delays + stat_function(fun=dgamma, 
    alpha=.5, args=list(shape=g_gg_alpha[i], scale=g_gg_beta[i]))

# How's it look:
pdf(file='gamma-model-st-lt-comparison.pdf', width=6, height=8)
gridExtra::grid.arrange(pl_g_st_delays, pl_g_lt_delays)
dev.off()

# That suggests the gamma model is not doing so good with longer tailed
# data... but it don't look so bad... how bad is it? Let's look at a q-q plot:
g_st_delays_tq <- vector(mode='numeric', length=length(st_delays_q))
for (i in seq_along(delay_quantiles)) {
  g_st_delays_tq[i] <- qgamma(p=delay_quantiles[i],
    shape=mean(g_g_alpha), scale=mean(g_g_beta))
}

g_lt_delays_tq <- vector(mode='numeric', length=length(lt_delays_q))
for (i in seq_along(delay_quantiles)) {
  g_lt_delays_tq[i] <- qgamma(p=delay_quantiles[i],
    shape=mean(g_gg_alpha), scale=mean(g_gg_beta))
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
if (!exists('generalized_gamma_on_gamma')) {
  generalized_gamma_on_gamma <- sampling(compiled_generalized_gamma, data=list(n_obs=N, n_parameters=3,
    t_start=rep(0, N), t_stop=st_delays, t_truncate=rep(Inf, N),
    kappa=rep(1,3), sigma_log=1/2, sigma_identity=1, sigma_logit=2/3),
    chains=5, iter=800, sample_file=paste0('generalized-gamma-on-gamma-out.'),
    diagnostic_file=paste0('generalized-gamma-on-gamma-diagnostic.')
  )
}


# Next we use the generalized gamma model on the heavier-tailed dataset:
if (!exists('generalized_gamma_on_generalized_gamma')) {
  generalized_gamma_on_generalized_gamma <- sampling(compiled_generalized_gamma, data=list(n_obs=N, n_parameters=3,
    t_start=rep(0, N), t_stop=lt_delays, t_truncate=rep(Inf, N),
    kappa=rep(1,3), sigma_log=1/2, sigma_identity=1, sigma_logit=2/3),
    chains=5, iter=800, sample_file=paste0('generalized-gamma-on-generalized-gamma-out.'),
    diagnostic_file=paste0('generalized-gamma-on-generalized-gamma-diagnostic.')
  )
}
# Next we use the generalized gamma model on the heavier-tailed dataset:
if (!exists('generalized_gamma_on_generalized_gamma')) {
  generalized_gamma_on_generalized_gamma <- sampling(compiled_generalized_gamma, data=list(n_obs=N, n_parameters=3,
    t_start=rep(0, N), t_stop=lt_delays, t_truncate=rep(Inf, N),
    kappa=rep(1,3), sigma_log=1/2, sigma_identity=1, sigma_logit=2/3),
    chains=5, iter=800, sample_file=paste0('generalized-gamma-on-generalized-gamma-out.'),
    diagnostic_file=paste0('generalized-gamma-on-generalized-gamma-diagnostic.')
  )
}

# Lets do a visual comparison of the two fits:
gg_g_alpha <- rstan::extract(generalized_gamma_on_gamma)[['alpha_0']] %>% exp
gg_g_beta <- rstan::extract(generalized_gamma_on_gamma)[['beta_0']] %>% exp
gg_g_nu <- rstan::extract(generalized_gamma_on_gamma)[['nu_0']] %>% exp
pl_gg_st_delays <- ggplot(data=data %>% filter(tail=='short'), 
  aes(x=delay)) + geom_histogram(aes(delay, ..density..), bins=100) +
  scale_x_continuous("delay", limits=c(0,8))
for (i in sample(x=1001:2000, size=10)) 
  pl_gg_st_delays <- pl_gg_st_delays + stat_function(fun=flexsurv::dgengamma.orig, alpha=.5, 
    args=list(shape=gg_g_nu[i], scale=gg_g_beta[i], k=gg_g_alpha[i]/gg_g_nu[i]))

gg_gg_alpha <- rstan::extract(generalized_gamma_on_generalized_gamma)[['alpha_0']] %>% exp
gg_gg_beta <- rstan::extract(generalized_gamma_on_generalized_gamma)[['beta_0']] %>% exp
gg_gg_nu <- rstan::extract(generalized_gamma_on_generalized_gamma)[['nu_0']] %>% exp
pl_gg_lt_delays <- ggplot(data=data %>% filter(tail=='long'), 
  aes(x=delay)) + geom_histogram(aes(delay, ..density..), bins=100) +
  scale_x_continuous("delay", limits=c(0,8))
for (i in sample(x=1001:2000, size=10)) 
  pl_gg_lt_delays <- pl_gg_lt_delays + stat_function(fun=flexsurv::dgengamma.orig, alpha=.5, 
    args=list(shape=gg_gg_nu[i], scale=gg_gg_beta[i], k=gg_gg_alpha[i]/gg_gg_nu[i]))

# How's it look:
pdf(file='generalized-gamma-model-st-lt-comparison.pdf', width=6, height=8)
gridExtra::grid.arrange(pl_gg_st_delays, pl_gg_lt_delays)
dev.off()

# That suggests the gamma model is not doing so good with longer tailed
# data... but it don't look so bad... how bad is it? Let's look at a q-q plot:
gg_st_delays_tq <- vector(mode='numeric', length=length(st_delays_q))
for (i in seq_along(delay_quantiles)) {
  gg_st_delays_tq[i] <- flexsurv::qgengamma.orig(p=delay_quantiles[i],
    shape=mean(gg_g_nu), scale=mean(gg_g_beta), k=mean(gg_g_alpha/gg_g_nu))
}

gg_lt_delays_tq <- vector(mode='numeric', length=length(lt_delays_q))
for (i in seq_along(delay_quantiles)) {
  gg_lt_delays_tq[i] <- flexsurv::qgengamma.orig(p=delay_quantiles[i],
    shape=mean(gg_gg_nu), scale=mean(gg_gg_beta), k=mean(gg_gg_alpha/gg_gg_nu))
}

gg_quantiles <- rbind(st_delays_q, gg_st_delays_tq, lt_delays_q, gg_lt_delays_tq) %>%
  data.frame(check.names=FALSE) %>% mutate(
    source=c(rep('short',2), rep('long',2)), 
    type=c('data','model','data','model')) %>%
  tidyr::gather(percentile, value, -source, -type) %>%
  tidyr::spread(type, value)

pl_gg_st_lt_qq_comparison <- ggplot(data=gg_quantiles, 
  aes(x=data, y=model, colour=percentile)
) + geom_point() +
    facet_wrap( ~ source, ncol=1) +
    geom_abline(slope=1, intercept=0) + 
    theme_minimal() + guides(colour=FALSE) +
    scale_x_continuous("delay (data quantiles)", limits=c(0,2.5)) +
    scale_y_continuous("delay (theoretical quantiles)", limits=c(0,2.5)) +
    coord_fixed()

# How's it look:
pdf(file='generalized-gamma-model-st-lt-qq-comparison.pdf', width=6, height=8)
print(pl_gg_st_lt_qq_comparison)
dev.off()

# Dump these to screen also:
gridExtra::grid.arrange(pl_gg_st_delays, pl_gg_lt_delays)
print(pl_gg_st_lt_qq_comparison)


# Four-way comarison by densities
pdf(file='model-st-lt-density-comparison.pdf', width=6, height=8)
gridExtra::grid.arrange(pl_g_st_delays, pl_g_lt_delays,
  pl_gg_st_delays, pl_gg_lt_delays, ncol=2)
dev.off()


all_quantiles <- rbind(
  g_quantiles %>% mutate(fit='gamma'),
  gg_quantiles %>% mutate(fit='generalized-gamma')
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




