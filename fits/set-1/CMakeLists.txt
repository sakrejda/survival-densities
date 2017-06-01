
add_custom_command(
  OUTPUT delay-data.rds
  COMMAND Rscript -e "source('${CMAKE_CURRENT_SOURCE_DIR}/data.R')"
  VERBATIM
)

add_custom_target(
  set-1-delay-data ALL DEPENDS delay-data.rds data.R
)


if (HAS_SCHEDULER)
  add_custom_command(
    OUTPUT model-binary-files.rds model-program-files.rds
    COMMAND Rscript ${CMAKE_CURRENT_SOURCE_DIR}/generate-submit-scripts.R
    DEPENDS delay-data.rds 
  )
else (HAS_SCHEDULER)
  add_custom_command(
    OUTPUT model-binary-files.rds model-program-files.rds
    COMMAND Rscript ${CMAKE_CURRENT_SOURCE_DIR}/generate-submit-scripts.R
    COMMAND Rscript ${CMAKE_CURRENT_SOURCE_DIR}/run.R job-files.txt ${CMAKE_CURRENT_SOURCE_DIR}
    DEPENDS delay-data.rds 
  )
endif (HAS_SCHEDULER)

add_custom_target(
  set-1-runs ALL DEPENDS model-binary-files.rds model-program-files.rds set-1-delay-data)

