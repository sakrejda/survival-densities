
add_custom_command(
  OUTPUT delay-data.rds
  COMMAND Rscript -e "source('${CMAKE_CURRENT_SOURCE_DIR}/good-data-feelings-data.R')"
  VERBATIM
)

add_custom_target(
  delay-data ALL DEPENDS delay-data.rds good-data-feelings-data.R
)

add_custom_command(
  OUTPUT data-quantiles.rds
  COMMAND Rscript -e "source('${CMAKE_CURRENT_SOURCE_DIR}/calculate-data-quantiles.R')"
  VERBATIM
)

add_custom_target(
  data-quantiles ALL DEPENDS data-quantiles.rds calculate-data-quantiles.R delay-data
)


add_custom_command(
  OUTPUT gamma-on-gamma-data.rdump gamma-on-gesgm-data.rdump generalized-gamma-on-gamma-data.rdump 
    generalized-gamma-on-gesgm-data.rdump gesgm-on-gamma-data.rdump gesgm-on-gesgm-data.rdump
  COMMAND Rscript -e "source('${CMAKE_CURRENT_SOURCE_DIR}/prep-model-data.R')"
  VERBATIM
)

add_custom_target(
  model-data-clean-contaminated ALL DEPENDS delay-data gamma-on-gamma-data.rdump gamma-on-gesgm-data.rdump generalized-gamma-on-gamma-data.rdump 
    generalized-gamma-on-gesgm-data.rdump gesgm-on-gamma-data.rdump gesgm-on-gesgm-data.rdump prep-model-data.R
)


add_custom_command(
  OUTPUT gesgm-on-gamma-inits.rdump gesgm-on-gesgm-inits.rdump
  COMMAND Rscript -e "source('${CMAKE_CURRENT_SOURCE_DIR}/prep-model-inits.R')"
  VERBATIM
)

add_custom_target(
  model-inits-clean-contaminated ALL DEPENDS gesgm-on-gamma-inits.rdump gesgm-on-gesgm-inits.rdump
)


if (IS_CLUSTER)
  set(FIND_ARGUMENTS ". -name \"\*-job-tag.pbs\" -exec ${SCHEDULER} < {} \\$<SEMICOLON>")
  separate_arguments(FIND_ARGUMENTS)
  add_custom_command(
    OUTPUT model-binary-files.rds model-program-files.rds
    COMMAND rm -f ${CMAKE_CURRENT_BINARY_DIR}/*-job-tag.pbs
    COMMAND Rscript ${CMAKE_CURRENT_SOURCE_DIR}/generate-submit-scripts.R
    COMMAND find ${FIND_ARGUMENTS}
    DEPENDS delay-data.rds stannis-build
  )
  add_custom_target(
    fat-tail-example-model-runs ALL DEPENDS model-binary-files.rds model-program-files.rds)
endif (NOT ${SCHEDULER} STREQUAL SCHEDULER-NOTFOUND)


