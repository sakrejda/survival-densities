library(simachine); library(magrittr); library(flexsurv)
sim_def <- yaml::yaml.load_file(input='parameters.yaml')

parameters <- fixed_parameter_points(sim_def$parameters)

saveRDS(object=parameters, file='N-alpha-beta-nu.rds')

