

add_custom_command(
  OUTPUT functions.stan
  COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/merge-stan-functions.sh ${CMAKE_CURRENT_SOURCE_DIR}
  VERBATIM)

add_custom_target(
  stan-functions ALL DEPENDS functions.stan)

add_custom_command(
  OUTPUT ${PROJECT_BINARY_DIR}/code/waitup/waitup/src/functions.hpp
  COMMAND make -C ${CMDSTAN_DIR} ${PROJECT_BINARY_DIR}/models/stan-lang/functions/functions
  COMMAND cp ${CMAKE_CURRENT_BINARY_DIR}/functions.hpp ${PROJECT_BINARY_DIR}/code/waitup/waitup/src/
  DEPENDS functions.stan
  VERBATIM)

add_custom_target(
  stan-functions-hpp ALL 
  DEPENDS ${PROJECT_BINARY_DIR}/code/waitup/waitup/src/functions.hpp waitup-dir stanc-binary)


