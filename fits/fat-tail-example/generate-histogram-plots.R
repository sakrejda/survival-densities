library(dplyr); library(ggplot2); library(magrittr); 

delay_data <- readRDS('delay-data.rds') %>% data.frame
model_densities <- readRDS('model-densities.rds')

pl_density <- list(
  clean = list(
    gamma = ggplot() + 
      geom_histogram(data=delay_data, aes(x=clean, y=..density..), bins=100) + 
      geom_line(data=model_densities$gamma$clean, aes(x=density, y=value, colour=factor(estimate))),
    gesgm = ggplot() + 
      geom_histogram(data=delay_data, aes(x=clean, y=..density..), bins=100) + 
      geom_line(data=model_densities$gesgm$clean, aes(x=density, y=value, colour=factor(estimate))),
    generalized_gamma = ggplot() + 
      geom_histogram(data=delay_data, aes(x=clean, y=..density..), bins=100) + 
      geom_line(data=model_densities$generalized_gamma$clean, aes(x=density, y=value, colour=factor(estimate)))
  ), 
  contaminated = list(
    generalized_gamma = ggplot() + 
      geom_histogram(data=delay_data, aes(x=contaminated, y=..density..), bins=100) + 
      geom_line(data=model_densities$generalized_gamma$contaminated, aes(x=density, y=value, colour=factor(estimate))),
    gamma = ggplot() + 
      geom_histogram(data=delay_data, aes(x=contaminated, y=..density..), bins=100) + 
      geom_line(data=model_densities$gamma$contaminated, aes(x=density, y=value, colour=factor(estimate))),
    gesgm = ggplot() + 
      geom_histogram(data=delay_data, aes(x=contaminated, y=..density..), bins=100) + 
      geom_line(data=model_densities$gesgm$contaminated, aes(x=density, y=value, colour=factor(estimate)))
  )
)

pl_density <- lapply(pl_density, function(x) lapply(x, function(x) {
  pl <- x + theme_minimal() + scale_x_continuous("delay") + 
    scale_y_continuous("density") 
}))

#grid.arrange(quantiles_plots$clean$gamma, quantiles_plots$clean$gesgm, quantiles_plots$contaminated$gamma, quantiles_plots$contaminated$gesgm, ncol=2)
delay_data_long <- delay_data %>% tidyr::gather(data, delay) %>% 
  rbind(.,.,.) %>% mutate(
    model=c(rep('gamma', 2*nrow(delay_data)),
            rep('generalized gamma', 2*nrow(delay_data)),
            rep('GESGM', 2*nrow(delay_data)))
  )

density_comparison <- list(
  model_densities$gamma$clean %>% mutate(data='clean', model='gamma'),
  model_densities$generalized_gamma$clean %>% mutate(data='clean', model='generalized gamma'),
  model_densities$gesgm$clean %>% mutate(data='clean', model='GESGM'),
  model_densities$gamma$contaminated %>% mutate(data='contaminated', model='gamma'),
  model_densities$generalized_gamma$contaminated %>% mutate(data='contaminated', model='generalized gamma'),
  model_densities$gesgm$contaminated %>% mutate(data='contaminated', model='GESGM')
) %>% do.call(what=rbind, args=.) %>% rename(model_d=value) #%>% 
#  left_join(y=delay_data %>% tidyr::gather(data, delay), by='data') 

pl_density_comparison <- ggplot() +
  geom_histogram(data=delay_data_long, aes(x=delay, y=..density..), bins=100) +
  geom_line(data=density_comparison, aes(x=density, y=model_d, colour=estimate), linetype=2) +
  facet_wrap( ~ data + model, ncol=3) + theme_minimal() +
  scale_x_continuous("delay", limits=c(0,20)) + 
  scale_y_continuous("density") 

pl_density_comparison_ur <- ggplot() +
  geom_histogram(data=delay_data_long, aes(x=delay, y=..density..), bins=100) +
  geom_line(data=density_comparison, aes(x=density, y=model_d, colour=estimate), linetype=2) +
  facet_wrap( ~ data + model, ncol=3) + theme_minimal() +
  scale_x_continuous("delay") + 
  scale_y_continuous("density") 

saveRDS(list(comparison=density_comparison, delays=delay_data_long), file='density-data.rds')

saveRDS(list(individual=pl_density, comparison=pl_density_comparison, 
             full_comparison=pl_density_comparison_ur), file='histogram-plots.rds')









