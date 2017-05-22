library(magrittr); library(digest)

parameters <- readRDS(file='N-alpha-beta-nu.rds')

data_targets <- list(
  `gamma-p1` = c('N', 'alpha', 'beta'),
  `weibull-p1` = c('N', 'alpha', 'beta'),
  `generalized-gamma-p1` = c('N', 'alpha', 'beta', 'nu')
)

target_parameters <- list()

for ( i in 1:length(data_targets)) {
  data_name <- names(data_targets)[i]
  duplicates <- lapply(X=parameters,
    FUN=function(x, names) x[names(x) %in% names],
    names=c(data_targets[[i]])) %>% sapply(digest) %>%
  duplicated

  target_parameters[[data_name]] <- lapply(X=parameters,
    FUN=function(x, names) x[names(x) %in% names],
    names=c(data_targets[[i]], 'point_digest'))
  target_parameters[[data_name]] <-
    target_parameters[[data_name]][!duplicates]

  saveRDS(object=target_parameters[[data_name]],
    file=paste0(data_name, '-data-target.rds'))
}
