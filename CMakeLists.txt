cmake_minimum_required (VERSION 3.2.2)
project (survival-densities VERSION 0.0.0.9000)

set(TOOLING_DIR ~/output-store/tooling)

set(CMDSTAN_DIR ${TOOLING_DIR}/code/stan-dev/cmdstan)
set(CMDSTAN_FIX_DIR ${TOOLING_DIR}/code/stan-dev/bugfix-cmdstan)
set(STAN_MODEL_DIR ${TOOLING_DIR}/models/stan-lang/full)
#set(ENV{CC} /share/pkg/gcc/5.1.0/bin/gcc)
#set(ENV{CXX} /share/pkg/gcc/5.1.0/bin/g++)
set(R_LIB_DIR ${TOOLING_DIR}/R/library)
set(ENV{R_LIBS_USER} ${R_LIB_DIR})

set(SIMULATION_DIR ${PROJECT_BINARY_DIR}/simulation)
set(FIT_DIR ${PROJECT_BINARY_DIR}/fit)

find_program(BSUB bsub)
if (NOT ${BSUB} STREQUAL BSUB-NOTFOUND)
  set(IS_BSUB TRUE)
  set(SCHEDULER "bsub")
else (NOT ${BSUB} STREQUAL BSUB-NOTFOUND)
  set(IS_BSUB FALSE)
endif(NOT ${BSUB} STREQUAL BSUB-NOTFOUND)

find_program(QSUB qsub)
if (NOT ${QSUB} STREQUAL QSUB-NOTFOUND)
  set(IS_QSUB TRUE)
  set(SCHEDULER "qsub")
else (NOT ${QSUB} STREQUAL QSUB-NOTFOUND)
  set(IS_QSUB FALSE)
endif (NOT ${QSUB} STREQUAL QSUB-NOTFOUND)

if (IS_BSUB OR IS_QSUB)
  set(HAS_SCHEDULER TRUE)
endif (IS_BSUB OR IS_QSUB)


#add_subdirectory(manuscripts)
add_subdirectory(simulation)
add_subdirectory(fits)
add_subdirectory(analysis)

