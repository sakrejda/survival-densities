
add_custom_command(
  OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/weibull-uniform-mixture.rds
  COMMAND cp ${CMAKE_CURRENT_SOURCE_DIR}/parameters.yaml ${CMAKE_CURRENT_BINARY_DIR}
  COMMAND "R_LIBS_USER=${R_LIB_DIR}" Rscript -e "source('${CMAKE_CURRENT_SOURCE_DIR}/sim.R')"
  DEPENDS parameters.yaml sim.R
  VERBATIM
)

add_custom_target(
  weibull-uniform-mixture-simulation ALL 
  DEPENDS ${CMAKE_CURRENT_BINARY_DIR}/weibull-uniform-mixture.rds
)



