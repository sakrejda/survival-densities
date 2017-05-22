library(magrittr); library(simachine); library(flexsurv)
source("../../../utilities.R")

sim_def <- yaml::yaml.load_file(input='parameters.yaml')

data_target_files <- dir(path='.', pattern='-data-target.rds')
data_target_names <- gsub(pattern='-data-target.rds', replacement='', x=data_target_files)

sim_code <- list(
  `gamma-p1` = quote({
    delay <- vector(mode='numeric', length=10^3)
    delay[1:ceiling(N)] <- rgamma(n=ceiling(N), shape=alpha, scale=beta)
  }),
  `weibull-p1` = quote({
    delay <- vector(mode='numeric', length=10^3)
    delay[1:ceiling(N)] <- rweibull(n=ceiling(N), shape=alpha, scale=beta)
  }),
  `generalized-gamma-p1` = quote({
    delay <- vector(mode='numeric', length=10^3)
    delay[1:ceiling(N)] <- rgengamma.orig(n=ceiling(N), shape=nu, scale=beta, k=alpha/nu)
  })
)

for ( name in data_target_names ) {
  parameters <- readRDS(file=paste0(name, '-data-target.rds'))
  parameter_names <- names(parameters[[1]])

  sim_output <- replicates(
    simulation=sim_code[[name]],
    parameters=parameters,
    control=list(
      size=sim_def$replicates$size,
      output=c('delay', parameter_names)
    )
  )

  for (i in 1:sim_def$replicates$size) {
    for (j in 1:length(parameters)) {
      data_set_name <- paste0(
        'data-target-', name, '--', 
        'theta-pt-', j, '--', 
        'data-set-', i)
      sim_result <- list(
        delay = sim_output$delay[,j,i],
        parameters = parameters[[j]]
      )
    
      saveRDS(object=sim_result, file=file.path('data', 
        paste0(data_set_name, '--simulation.rds'))) 
      
      N <- sim_result$parameters$N
      sim_data <- stan_data(t_start=rep(0, N), t_stop=sim_result$delay[1:N], 
        kappa=list(alpha=1, beta=1, nu=1, k=1, mu=1, sigma=1))
      rstan::stan_rdump(list=names(sim_data), file=file.path('data', 
        paste0(data_set_name, '--simulation.rdump')), envir=list2env(sim_data))

    }
  }
}



