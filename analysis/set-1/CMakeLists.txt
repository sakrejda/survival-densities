
add_custom_command(
  OUTPUT run_metadata.rds
  COMMAND "R_LIBS_USER=${R_LIB_DIR}" Rscript -e "source('${CMAKE_CURRENT_SOURCE_DIR}/run_metadata.R')"
  VERBATIM
)

add_custom_target(
  set-1-completion ALL DEPENDS set-1-runs run_metadata.rds
)

