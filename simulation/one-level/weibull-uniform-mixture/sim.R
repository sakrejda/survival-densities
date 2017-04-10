library(simachine);
sim_def <- yaml::yaml.load_file(input='parameters.yaml')

parameters <- random_parameter_trees(sim_def$parameters,
  sim_def$replicates$draws_per_parameter)

sim_code <- quote({
  z <- runif(n=ceiling(N))
  ru <- runif(n=ceiling(N), min=0, max=delay_max)
  rw <- rweibull(n=ceiling(N), shape=alpha, scale=beta)
  delay <- rep(NA, 10^3)
  delay[1:ceiling(N)] <- ifelse(z < q, ru, rw)
})

sim_output <- replicates(
  simulation=sim_code,
  parameters=parameters,
  control=list(
    size=sim_def$replicates$size,
    output=c('delay','N','alpha','beta','q','delay_max')
  )
)

sim_results <- list(
  delay=sim_output$delay, 
  parameters=with(data=sim_output, expr=abind::abind(
    N=N, alpha=alpha, beta=beta, q=q, delay_max=delay_max, along=1))
)

saveRDS(object=sim_results, file='weibull-uniform-mixture.rds')

# 4) (later fit a model to this density)
# 5) recover the parameters
# 6) compare to the density


## Add installation CMakeLists for: yaml, simachine, dplyr, ggplot,
## abind
## etc...



