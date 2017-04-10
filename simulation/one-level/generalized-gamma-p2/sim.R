library(simachine);
sim_def <- yaml::yaml.load_file(input='parameters.yaml')

parameters <- random_parameter_trees(sim_def$parameters,
  sim_def$replicates$draws_per_parameter)

sim_code <- quote({
  delay <- rep(NA, 10^3)
  delay[1:ceiling(N)] <- flexsurv::rgengamma(n=ceiling(N), mu=mu, sigma=sigma, Q=1/sqrt(k)) 
})

sim_output <- replicates(
  simulation=sim_code,
  parameters=parameters,
  control=list(
    size=sim_def$replicates$size,
    output=c('delay','N','k','mu','sigma')
  )
)

sim_results <- list(
  delay=sim_output$delay,
  parameters=with(data=sim_output, abind::abind(
    N=N, k=k, mu=mu, sigma=sigma, along=1))
)

saveRDS(object=sim_results, file='generalized-gamma-p2.rds')

# 4) (later fit a model to this density)
# 5) recover the parameters
# 6) compare to the density


## Add installation CMakeLists for: yaml, simachine, dplyr, ggplot,
## abind
## etc...



