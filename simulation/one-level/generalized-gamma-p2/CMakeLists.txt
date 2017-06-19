
add_custom_command(
  OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/generalized-gamma-p2.rds
  COMMAND cp ${CMAKE_CURRENT_SOURCE_DIR}/parameters.yaml ${CMAKE_CURRENT_BINARY_DIR}
  COMMAND "R_LIBS_USER=${R_LIB_DIR}" Rscript -e "source('${CMAKE_CURRENT_SOURCE_DIR}/sim.R')"
  DEPENDS parameters.yaml sim.R
  VERBATIM
)

add_custom_target(
  generalized-gamma-p2-simulation ALL 
  DEPENDS ${CMAKE_CURRENT_BINARY_DIR}/generalized-gamma-p2.rds
)



