library(magrittr)

output_path <- '../../fits/set-1'

run_data <- readRDS(file=file.path(output_path, 'run-data.rds'))
n_runs <- nrow(run_data)

run_data <- dplyr::mutate(run_data, 
  completed_iterations = 0,
  completed_warmup = FALSE,
  sampling_stepsize = NULL,
  sampling_divergence_count = NULL,
  sampling_n_leapfrog_mean = NULL,
  sampling_n_leapfrog_sd = NULL,
  single_chain_r_hat_mean = NULL,
  single_chain_r_hat_sd = NULL,
  total_time = NA
)

all_estimates <- NULL

for (i in 1:n_runs) {
  output_file <- file.path(output_path, run_data[i, 'output_file'])
  if (!file.exists(output_file))
    next
  messed_up <- try(expr = {
    lines <- readLines(output_file)
    samples <- stannis::read_stan_data(file=output_file)
    meta <- stannis::read_stan_metadata(output_file)
  }, silent=TRUE) 
  if (class(messed_up) == "try-error")
    next
  run_data[i, 'completed_iterations'] <- nrow(samples)
  run_data[i, 'completed_warmup'] <- nrow(samples) > meta[['num_warmup']]
  run_data[i, 'completed_sampling'] <- 
    nrow(samples) == (meta[['num_warmup']] + meta[['num_draws']] - 1)
  if ( run_data[i, 'completed_sampling'] ) {
    sampling <- (meta[['num_warmup']]+1):nrow(samples)
    run_data[i, 'sampling_stepsize'] <- unique(samples[sampling, 'stepsize__'])
    run_data[i, 'divergence_count'] <- sum(samples[['divergent__']])
    run_data[i, 'sampling_n_leapfrog_mean'] <- mean(samples[['n_leapfrog__']])
    run_data[i, 'sampling_n_leapfrog_sd'] <- sd(samples[['n_leapfrog__']])
    if (nrow(samples) > 500) {
      r_hat_samples <- samples[sampling, c(1, 8:ncol(samples))]
      run_data[i, 'single_chain_r_hat_mean'] <- sapply(r_hat_samples,
        rstan:::split_rhat_rfun) %>% mean(na.rm=TRUE)
      run_data[i, 'single_chain_r_hat_sd'] <- sapply(r_hat_samples, 
        rstan:::split_rhat_rfun) %>% sd(na.rm=TRUE)
    } else {
      run_data[i, 'single_chain_r_hat_mean'] <- 2
      run_data[i, 'single_chain_r_hat_sd'] <- 2
    }
    if (nrow(samples) == 1500) {
      total_time <- stannis::get_total_time(lines)
      if (!is.na(total_time) && total_time < 30*60) 
        run_data[i, 'total_time'] <- total_time
      else
        run_data[i, 'total_time'] <- NA  
    }
    if (grepl(pattern='sample', x=meta[['method']])) {
      estimates <- rstan::read_stan_csv(output_file) %>% 
        as.matrix %>% 
        apply(2, quantile, probs=c(0.05, 0.25, 0.5, 0.75, 0.95)) %>% 
        t %>% data.frame(check.names=FALSE) %>% 
        dplyr::mutate(model=run_data[i, 'model'], data=run_data[i, 'data'], 
          binary_type=run_data[i, 'binary_type'], parameter=rownames(.),
          chain=run_data[i, 'chain'])
      if (is.null(all_estimates)) 
        all_estimates <- estimates
      else 
        all_estimates <- rbind(all_estimates, estimates)
    }
  }
}

saveRDS(run_data, file='run_metadata.rds')

warmup_completion <- run_data %>% 
  dplyr::group_by(model, data_type) %>% dplyr::summarise(
    percent_completed=mean(completed_warmup)*100
  )

sampling_completion <- run_data %>% 
  dplyr::group_by(model, data_type) %>% dplyr::summarise(
    percent_completed=mean(completed_sampling)*100
  )

timing <- run_data %>% 
  dplyr::group_by(model, data_type) %>% dplyr::summarise(
    mean_time=mean(total_time, na.rm=TRUE),
    sd_time=sd(total_time, na.rm=TRUE)
  )

run_type <- run_data %>% dplyr::filter(completed_sampling) %>% 
  dplyr::select(model, data_type output_file) %>%
  dplyr::group_by(model, data_type) %>% 
  dplyr::summarise(output_file=paste(file.path(output_path, output_file), collapse=' ')) %>%
  lapply(as.character) %>% data.frame(stringsAsFactors=FALSE)

for (i in 1:nrow(run_type)) {
  binary <- '../../code/stan-dev/cmdstan/bin/stansummary'
  files <- run_type[i, 'output_file']
  out <- paste0(run_type[i, 'model'], "-on-", run_type[i, 'data'],
    "-with-", run_type[i, 'binary_type'], "-stan-summary.csv")
  if (file.exists(out)) file.remove(out)
  cmd <- paste(binary, files, paste0("--csv_file=", out))
  system(cmd)
  if (!file.exists(out)) next
  run_type_summary <- read.csv(out, check.names=FALSE, comment.char="#", 
    stringsAsFactors=FALSE)
  non_sampler_parameters <- 
    run_type_summary[['name']][!grepl(pattern='__$', x=run_type_summary[['name']])]
  run_type[i, 'multi_chain_R_hat_mean'] <- 
    dplyr::filter(run_type_summary, name %in% c('lp__', non_sampler_parameters)) %>%
    dplyr::select(R_hat) %>% unlist %>% mean
  run_type[i, 'multi_chain_R_hat_sd'] <- 
    dplyr::filter(run_type_summary, name %in% c('lp__', non_sampler_parameters)) %>%
    dplyr::select(R_hat) %>% unlist %>% sd
}


pl_scaled_parameters <- ggplot(
  data=all_estimates %>% dplyr::filter(
    parameter %in% c('g_alpha', 'g_beta', 'g_delta') & 
    (model != 'gamma-exp-sum-gamma-mix-p2' | data!='gesgm' | binary_type!='fix' | chain != 5)), 
  aes(xmin=`25%`, x=`50%`, xmax=`95%`, y=model, colour=binary_type)
) + geom_point() + 
    geom_errorbarh() + 
    facet_grid(data_type ~ parameter, scales='free_x')

pl_internal_parameters <- ggplot(
  data=all_estimates %>% dplyr::filter(
    parameter %in% c('mu_0', 'sigma_0', 'alpha_0', 'beta_0', 'delta_0')), 
  aes(xmin=`25%`, x=`50%`, xmax=`95%`, y=model, colour=binary_type)
) + geom_point() + 
    geom_errorbarh() + 
    facet_grid(data_type ~ parameter, scales='free_x')


