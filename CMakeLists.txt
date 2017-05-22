cmake_minimum_required (VERSION 3.2.2)
project (delay-analysis-method VERSION 0.0.0.9000)

set(CMDSTAN_DIR ${PROJECT_BINARY_DIR}/code/stan-dev/cmdstan)
set(CMDSTAN_FIX_DIR ${PROJECT_BINARY_DIR}/code/stan-dev/bugfix-cmdstan)
set(SIMULATION_DIR ${PROJECT_BINARY_DIR}/simulation)
set(FIT_DIR ${PROJECT_BINARY_DIR}/fit)
set(STAN_MODEL_DIR ${PROJECT_BINARY_DIR}/models/stan-lang/full)
set(STAN_MODEL_FIX_DIR ${PROJECT_BINARY_DIR}/models/stan-lang/full)
#set(ENV{CC} /share/pkg/gcc/5.1.0/bin/gcc)
#set(ENV{CXX} /share/pkg/gcc/5.1.0/bin/g++)
set(R_LIB_DIR ${PROJECT_BINARY_DIR}/R/library)
set(ENV{R_LIBS_USER} ${R_LIB_DIR})

add_custom_command(
  OUTPUT ${PROJECT_BINARY_DIR}/utilities.R
  COMMAND cp ${PROJECT_SOURCE_DIR}/utilities.R ${PROJECT_BINARY_DIR}
  COMMAND touch ${PROJECT_BINARY_DIR}/utilities.done
)

add_custom_target(
  utilities-code ALL DEPENDS ${PROJECT_BINARY_DIR}/utilities.R
)

add_custom_command(
  OUTPUT R
  COMMAND mkdir -p ${R_LIB_DIR}
  COMMAND Rscript ${PROJECT_SOURCE_DIR}/packages.R ${PROJECT_SOURCE_DIR} noupdate
  VERBATIM
  DEPENDS packages.R R-libraries.txt simachine-build
)

add_custom_target(
  r-environment ALL DEPENDS R
)

add_subdirectory(code)
add_subdirectory(models)
#add_subdirectory(manuscripts)
add_subdirectory(simulation)
add_subdirectory(fits)

