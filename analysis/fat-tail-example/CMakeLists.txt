
add_custom_command(
  OUTPUT run_metadata.rds
  COMMAND "R_LIBS_USER=${R_LIB_DIR}" Rscript -e "source('${CMAKE_CURRENT_SOURCE_DIR}/run_metadata.R')"
  VERBATIM
)

add_custom_target(
  fat-tail-completion ALL DEPENDS fat-tail-runs run_metadata.rds
)

