
add_custom_command(
  OUTPUT N-alpha-beta-nu.rds
  COMMAND cp ${CMAKE_CURRENT_SOURCE_DIR}/parameters.yaml ${CMAKE_CURRENT_BINARY_DIR}
  COMMAND Rscript -e "source('${CMAKE_CURRENT_SOURCE_DIR}/make-parameter-points.R')"
  VERBATIM
)

add_custom_target(
  set-1-parameter-points ALL DEPENDS N-alpha-beta-nu.rds simachine-build
)

add_custom_command(
  OUTPUT gamma-p1-data-target.rds weibull-p1-data-target.rds generalized-gamma-p1-data-target.rds
  COMMAND Rscript -e "source('${CMAKE_CURRENT_SOURCE_DIR}/make-data-targets.R')"
  VERBATIM
)

add_custom_target(
  set-1-data-targets ALL DEPENDS gamma-p1-data-target.rds 
    weibull-p1-data-target.rds generalized-gamma-p1-data-target.rds
)

add_custom_command(
  OUTPUT data 
  COMMAND mkdir -p ${CMAKE_CURRENT_BINARY_DIR}/data
  COMMAND Rscript -e "source('${CMAKE_CURRENT_SOURCE_DIR}/make-data-sets.R')"
  VERBATIM
)

add_custom_target(
  set-1-data-sets ALL DEPENDS data 
)


