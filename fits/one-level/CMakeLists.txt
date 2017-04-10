add_subdirectory(set-1)

### There should be TWO targets here, one pointing to the random-targets
### data and one pointing to the set-1 data directory.  Both use 
### fit.R to run.

#set (SET-1-FITS  
#  set-1/gamma-p1--gamma-p1.rds set-1/gamma-p1--weibull-p1.rds
#  set-1/gamma-p1--generalized-gamma-p1.rds 
#  set-1/weibull-p1--gamma-p1.rds set-1/weibull-p1--weibull-p1.rds
#  set-1/weibull-p1--generalized-gamma-p1.rds 
#  set-1/generalized-gamma-p1--gamma-p1.rds set-1/generalized-gamma-p1--weibull-p1.rds
#  set-1/generalized-gamma-p1--generalized-gamma-p1.rds 
#)
#
#add_custom_command(
#  OUTPUT ${SET-1-FITS}
#  COMMAND mkdir -p ${CMAKE_CURRENT_BINARY_DIR}/set-1
#  COMMAND Rscript ${CMAKE_CURRENT_SOURCE_DIR}/fit.R ${PROJECT_BINARY_DIR}/models/stan-lang/full ${SIMULATION_DIR}/one-level/set-1 ${CMAKE_CURRENT_BINARY_DIR}/set-1
#  DEPENDS gamma-p1 weibull-p1
#    weibull-uniform-mixture-p1 generalized-gamma-p1 generalized-gamma-p2
#  VERBATIM
#)
#
#add_custom_target(set-1-output ALL DEPENDS set-1-simulations ${SET-1-FITS})

