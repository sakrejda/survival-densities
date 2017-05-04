library(dplyr); library(ggplot2); library(magrittr); 

data_quantiles <- readRDS('data-quantiles.rds')
model_quantiles <- readRDS('model-quantiles.rds')

quantiles_plots <- list(
  clean = list(
    gamma = ggplot() + 
      geom_line(data=data_quantiles$data$clean, aes(x=quantile, y=value)) + 
      geom_line(data=model_quantiles$gamma$clean, aes(x=quantile, y=value, colour=factor(estimate))),
    gesgm = ggplot() + 
      geom_line(data=data_quantiles$data$clean, aes(x=quantile, y=value)) + 
      geom_line(data=model_quantiles$gesgm$clean, aes(x=quantile, y=value, colour=factor(estimate))),
    generalized_gamma = ggplot() + 
      geom_line(data=data_quantiles$data$clean, aes(x=quantile, y=value)) + 
      geom_line(data=model_quantiles$generalized_gamma$clean, aes(x=quantile, y=value, colour=factor(estimate)))
  ), 
  contaminated = list(
    generalized_gamma = ggplot() + 
      geom_line(data=data_quantiles$data$contaminated, aes(x=quantile, y=value)) + 
      geom_line(data=model_quantiles$generalized_gamma$contaminated, aes(x=quantile, y=value, colour=factor(estimate))),
    gamma = ggplot() + 
      geom_line(data=data_quantiles$data$contaminated, aes(x=quantile, y=value)) + 
      geom_line(data=model_quantiles$gamma$contaminated, aes(x=quantile, y=value, colour=factor(estimate))),
    gesgm = ggplot() + 
      geom_line(data=data_quantiles$data$contaminated, aes(x=quantile, y=value)) + 
      geom_line(data=model_quantiles$gesgm$contaminated, aes(x=quantile, y=value, colour=factor(estimate)))
  )
)

quantiles_plots <- lapply(quantiles_plots, function(x) lapply(x, function(x) {
  pl <- x + theme_minimal() + scale_x_continuous("probability") + 
    scale_y_continuous("quantile") 
}))

#grid.arrange(quantiles_plots$clean$gamma, quantiles_plots$clean$gesgm, quantiles_plots$contaminated$gamma, quantiles_plots$contaminated$gesgm, ncol=2)


quantile_comparison <- list(
  model_quantiles$gamma$clean %>% mutate(data='clean', model='gamma'),
  model_quantiles$generalized_gamma$clean %>% mutate(data='clean', model='generalized gamma'),
  model_quantiles$gesgm$clean %>% mutate(data='clean', model='GESGM'),
  model_quantiles$gamma$contaminated %>% mutate(data='contaminated', model='gamma'),
  model_quantiles$generalized_gamma$contaminated %>% mutate(data='contaminated', model='generalized gamma'),
  model_quantiles$gesgm$contaminated %>% mutate(data='contaminated', model='GESGM')
) %>% do.call(what=rbind, args=.) %>% rename(model_q=value) %>%  
  mutate(quantile=as.character(quantile)) %>% left_join(
    y=data_quantiles$data %>% rename(data_q=value) %>% mutate(quantile=as.character(quantile)), 
    by=c('data', 'quantile'), copy=TRUE) %>% mutate(quantile=as.numeric(quantile))

pl_quantile_comparison <- ggplot(data=quantile_comparison) +
  geom_line(aes(x=quantile, y=data_q), linetype=1) + 
  geom_line(aes(x=quantile, y=model_q, colour=estimate), linetype=2) +
  facet_wrap( ~ data + model, ncol=3) + theme_minimal() +
  scale_x_continuous("probability") + 
  scale_y_continuous("quantile") 


saveRDS(list(comparison=quantile_comparison), file='qq-data.rds')

saveRDS(list(individual=quantiles_plots, comparison=pl_quantile_comparison), file='qq-plots.rds')









