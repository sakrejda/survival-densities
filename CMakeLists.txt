cmake_minimum_required (VERSION 3.5)
project (delay-analysis-method VERSION 0.0.0.9000)

set(CMDSTAN_DIR ${PROJECT_BINARY_DIR}/code/stan-dev/cmdstan)
set(CMDSTAN_FIX_DIR ${PROJECT_BINARY_DIR}/code/stan-dev/bugfix-cmdstan)
set(SIMULATION_DIR ${PROJECT_BINARY_DIR}/simulation)
set(FIT_DIR ${PROJECT_BINARY_DIR}/fit)
set(STAN_MODEL_DIR ${PROJECT_BINARY_DIR}/models/stan-lang/full)
set(R_LIB_DIR ${PROJECT_BINARY_DIR}/packrat/*/*)

add_custom_command(
  OUTPUT ${PROJECT_BINARY_DIR}/utilities.R
  COMMAND cp ${PROJECT_SOURCE_DIR}/utilities.R ${PROJECT_BINARY_DIR}
  COMMAND touch ${PROJECT_BINARY_DIR}/utilities.done
)

add_custom_target(
  utilities-code ALL DEPENDS ${PROJECT_BINARY_DIR}/utilities.R
)

add_custom_command(
  OUTPUT packrat
  COMMAND export R_LIBS_USER=${R_LIB_DIR}
  COMMAND Rscript ${PROJECT_SOURCE_DIR}/packrat.R ${PROJECT_SOURCE_DIR} noupdate
  VERBATIM
  DEPENDS packrat.R R-libraries.txt simachine-build
)

add_custom_target(
  r-environment ALL DEPENDS packrat
)

add_subdirectory(code)
add_subdirectory(models)
add_subdirectory(manuscripts)
add_subdirectory(simulation)
add_subdirectory(fits)

