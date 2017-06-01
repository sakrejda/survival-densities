
add_custom_command(
  OUTPUT delay-data.rds
  COMMAND Rscript -e "source('${CMAKE_CURRENT_SOURCE_DIR}/good-data-feelings-data.R')"
  VERBATIM
)

add_custom_target(
  delay-data ALL DEPENDS delay-data.rds good-data-feelings-data.R
)


if (IS_CLUSTER)
  set(FIND_ARGUMENTS ". -name \"\*-job-tag.pbs\" -exec ${SCHEDULER} < {} \\$<SEMICOLON>")
  separate_arguments(FIND_ARGUMENTS)
  add_custom_command(
    OUTPUT model-binary-files.rds model-program-files.rds
    COMMAND rm -f ${CMAKE_CURRENT_BINARY_DIR}/*-job-tag.pbs
    COMMAND Rscript ${CMAKE_CURRENT_SOURCE_DIR}/generate-submit-scripts.R
    COMMAND find ${FIND_ARGUMENTS}
    DEPENDS delay-data.rds 
  )
  add_custom_target(
    fat-tail-example-model-runs ALL DEPENDS model-binary-files.rds model-program-files.rds)
endif (IS_CLUSTER)


