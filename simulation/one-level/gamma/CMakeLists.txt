
add_custom_command(
  OUTPUT gamma.rds
  COMMAND cp ${CMAKE_CURRENT_SOURCE_DIR}/parameters.yaml ${CMAKE_CURRENT_BINARY_DIR}
  COMMAND Rscript -e "source('${CMAKE_CURRENT_SOURCE_DIR}/sim.R')"
  VERBATIM
)

add_custom_target(
  gamma-simulation ALL DEPENDS gamma.rds
)



