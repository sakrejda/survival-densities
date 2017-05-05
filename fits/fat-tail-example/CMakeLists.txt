
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
  data-quantiles ALL DEPENDS data-quantiles.rds calculate-data-quantiles.R
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

find_program(SCHEDULER bsub)

if (NOT ${SCHEDULER} STREQUAL SCHEDULER-NOTFOUND)
  set(FIND_ARGUMENTS ". -name \"\*-job.pbs\" -exec bsub < {} \\$<SEMICOLON>")
else (NOT ${SCHEDULER} STREQUAL SCHEDULER-NOTFOUND)
  set(FIND_ARGUMENTS ". -name \"\*-job.pbs\" -exec Rscript run.R {} \\$<SEMICOLON>")
endif (NOT ${SCHEDULER} STREQUAL SCHEDULER-NOTFOUND)

