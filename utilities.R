

drop_extension <- function(x) gsub('\\.[^\\.]*', x, replacement='')

stan_data <- function(
  t_start,
  t_stop,
  t_truncate = NULL,
  kappa = NULL,
  prior_sd = list(log=1/2, logit=3/2)
) {
  o <- list()
  if (is.null(t_truncate))
    t_truncate <- rep(Inf, length(t_stop))
  if (!all(length(t_start) == length(t_stop) &&
           length(t_start) == length(t_truncate)))
    stop("t_start, t_stop, and t_truncate must be the same length.")
  o[['t_start']] <- t_start
  o[['t_stop']] <- t_stop
  o[['t_truncate']] <- t_truncate
  o[['n_obs']] <- length(t_start)
  o[['n_parameters']] <- length(kappa)
  if (!is.null(kappa)) {
    o[['kappa']] <- vector(mode='numeric', length=length(kappa))
    for (i in seq_along(kappa))
      o[['kappa']][i] <- kappa[[i]]
  } 
  for (name in names(prior_sd))
    o[[paste('sigma', name, sep='_')]] <- prior_sd[[name]]
  return(o)
}


fit <- function(stan_path, data_path, output_path) {
  options(mc.cores=6)
  library(rstan)
  model_name <- basename(stan_path) %>% strsplit('\\.') %>% `[[`(1) %>% `[`(1)
  data_name <- basename(data_path) %>% strsplit('\\.') %>% `[[`(1) %>% `[`(1)
  output_name <- basename(output_path) %>% strsplit('\\.') %>% `[[`(1) %>% `[`(1)
  
  sim <- readRDS(file=data_path)
  data <- sim$delay
  theta <- sim$parameters
  n_data_replicates <- dim(data)[3]
  n_parameter_replicates <- dim(theta)[2]
  parameter_names <- dimnames(theta)[[1]]
  
  compiled_model <- stan_model(file=stan_path, model_name=model_name,
    save_dso=TRUE, obfuscate_model_name=FALSE)
  
  for(i in 1:n_parameter_replicates) {
    N <- ceiling(theta['N',i,1]) ## Same for all data replicates.
    parameter_replicate_name <- paste("theta-number", i, sep='-')
    for(j in 1:n_data_replicates) {
      t_stop <- data[,i,j,drop=TRUE][1:N]
      t_start <- rep(0, length(t_stop))
      formated_data <- stan_data(t_start=t_start, t_stop=t_stop,
        kappa=list(alpha=1, beta=1, nu=1, k=1, mu=1, sigma=1))
      data_replicate_name <- paste("data-number", j, sep='-')
      fit_name <- paste(model_name, data_name, 
        parameter_replicate_name, data_replicate_name, sep='_') 
      results <- sampling(compiled_model, data=formated_data, chains=5, iter=800,
        sample_file=paste0(fit_name,'.out'), 
        diagnostic_file=paste0(fit_name,'.diagnostic')
      )
      saveRDS(object=results, file=output_path)
    }
  }
}





