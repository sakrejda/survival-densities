
add_custom_command(
  OUTPUT completion.R
  COMMAND Rscript -e "source('${CMAKE_CURRENT_SOURCE_DIR}/good-data-feelings-data.R')"
  VERBATIM
)

add_custom_target(
  fat-tail-example-completion ALL DEPENDS fat-tail-example-model-runs completion.R
)

