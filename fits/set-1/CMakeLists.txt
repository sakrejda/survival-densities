
add_custom_command(
  OUTPUT delay-data.rds
  COMMAND Rscript -e "source('${CMAKE_CURRENT_SOURCE_DIR}/data.R')"
  VERBATIM
)

add_custom_target(
  set-1-delay-data ALL DEPENDS delay-data.rds data.R
)

find_program(SCHEDULER bsub)

if (NOT ${SCHEDULER} STREQUAL SCHEDULER-NOTFOUND)
  set(FIND_ARGUMENTS ". -name \"\*-job-tag.pbs\" -exec bsub < {} \\$<SEMICOLON>")
  separate_arguments(FIND_ARGUMENTS)
  add_custom_command(
    OUTPUT model-binary-files.rds model-program-files.rds
    COMMAND rm -f ${CMAKE_CURRENT_BINARY_DIR}/*-job-tag.pbs
    COMMAND Rscript ${CMAKE_CURRENT_SOURCE_DIR}/generate-submit-scripts.R
    COMMAND find ${FIND_ARGUMENTS}
    DEPENDS delay-data.rds stannis-build
  )
else (NOT ${SCHEDULER} STREQUAL SCHEDULER-NOTFOUND)
  set(FIND_ARGUMENTS_1 ". -name \"\*-job-tag.pbs\" -exec chmod 755 {} \\$<SEMICOLON>")
  set(FIND_ARGUMENTS_2 ". -name \"\*-job-tag.pbs\" > job-files.txt")
  separate_arguments(FIND_ARGUMENTS_1)
  separate_arguments(FIND_ARGUMENTS_2)
  add_custom_command(
    OUTPUT model-binary-files.rds model-program-files.rds
    COMMAND rm -f ${CMAKE_CURRENT_BINARY_DIR}/*-job-tag.pbs
    COMMAND Rscript ${CMAKE_CURRENT_SOURCE_DIR}/generate-submit-scripts.R
    COMMAND find ${FIND_ARGUMENTS_1}
    COMMAND find ${FIND_ARGUMENTS_2}
    COMMAND Rscript ${CMAKE_CURRENT_SOURCE_DIR}/run.R job-files.txt ${CMAKE_CURRENT_SOURCE_DIR}
    DEPENDS delay-data.rds stannis-build
  )
endif (NOT ${SCHEDULER} STREQUAL SCHEDULER-NOTFOUND)

add_custom_target(
  set-1-runs ALL DEPENDS model-binary-files.rds model-program-files.rds set-1-delay-data)

