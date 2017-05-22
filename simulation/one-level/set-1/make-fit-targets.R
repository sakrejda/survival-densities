library(magrittr)

parameters <- readRDS(file='N-alpha-beta-nu.rds')

fit_targets <- list(
  `gamma-p1` = c('N', 'alpha', 'beta'),
  `weibull-p1` = c('N', 'alpha', 'beta'),
  `generalized-gamma-p1` = c('N', 'alpha', 'beta', 'nu'),
  `generalized-gamma-p2` = c('N', 'k', 'mu', 'sigma')
)

data_target_files <- dir(path='.', pattern='-data-targets.rds')

for ( i in 1:length(data_targets)) {
  fit_target_name <- names(fit_targets)[i]
  

}
