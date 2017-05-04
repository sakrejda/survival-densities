library(magrittr)

my_rdump <- function(x, file) {
  rstan::stan_rdump(names(x), file, envir=list2env(x))
}

gamma_on_gamma_inits <- list(
) %>% my_rdump(file='gamma-on-gamma-inits.rdump')
  
gamma_on_gesgm_inits <- list(
) %>% my_rdump(file='gamma-on-gesgm-inits.rdump')

generalized_gamma_on_gamma_inits <- list(
) %>% my_rdump(file='generalized-gamma-on-gamma-inits.rdump')

generalized_gamma_on_gesgm_inits <- list(
) %>% my_rdump(file='generalized-gamma-on-gesgm-inits.rdump')

gesgm_on_gamma_inits <- list(
  mu_0 = sqrt(3), sigma_0=sqrt(sqrt(3)), q_0=-5, delta_0=0
) %>% my_rdump(file='gesgm-on-gamma-inits.rdump')

gesgm_on_gesgm_inits <- list(
  mu_0 =log(sqrt(3)*sqrt(3)), sigma_0=log(sqrt(3)), q_0=-1, delta_0=log(5/1)-3
) %>% my_rdump(file='gesgm-on-gesgm-inits.rdump')


