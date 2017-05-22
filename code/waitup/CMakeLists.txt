
set(ENV{R_LIBS_USER} ${R_LIB_DIR})

file (GLOB RPKG waitup/*)

add_custom_command(
  OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/waitup
  COMMAND rm -rf ${CMAKE_CURRENT_BINARY_DIR}/waitup
  COMMAND export R_LIBS_USER=${R_LIB_DIR}
  COMMAND cp ${CMAKE_CURRENT_SOURCE_DIR}/.Renviron ${CMAKE_CURRENT_BINARY_DIR}/.Renviron
  COMMAND cp -R ${CMAKE_CURRENT_SOURCE_DIR}/waitup ${CMAKE_CURRENT_BINARY_DIR}
  COMMAND mkdir -p ${CMAKE_CURRENT_BINARY_DIR}/waitup/inst/include
  COMMAND cp -R ${PROJECT_BINARY_DIR}/code/stan-dev/stan ${CMAKE_CURRENT_BINARY_DIR}/waitup/inst/include/
  DEPENDS ${RPKG} 
  VERBATIM
)

add_custom_target(
  waitup-dir ALL DEPENDS ${CMAKE_CURRENT_BINARY_DIR}/waitup stan-repo)

add_custom_command(
  OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/waitup_${PROJECT_VERSION}.tar.gz
  COMMAND export R_LIBS_USER=${R_LIB_DIR}
  COMMAND Rscript -e "Rcpp:::compileAttributes(pkgdir='waitup')"
  COMMAND R_LIBS_USER=${R_LIB_DIR} R CMD build ${CMAKE_CURRENT_BINARY_DIR}/waitup
  COMMAND R CMD INSTALL ${CMAKE_CURRENT_BINARY_DIR}/waitup_${PROJECT_VERSION}.tar.gz --library=${R_LIB_DIR}
  COMMAND Rscript -e ".libPaths(\"${R_LIB_DIR}\"); devtools::document('waitup')"
  COMMAND R_LIBS_USER=${R_LIB_DIR} R CMD build ${CMAKE_CURRENT_BINARY_DIR}/waitup
  COMMAND R CMD INSTALL ${CMAKE_CURRENT_BINARY_DIR}/waitup_${PROJECT_VERSION}.tar.gz --library=${R_LIB_DIR}
  DEPENDS waitup-dir stan-functions-hpp 
  VERBATIM)

add_custom_target(
  waitup-pkg ALL DEPENDS waitup-dir stan-functions-hpp 
    ${CMAKE_CURRENT_BINARY_DIR}/waitup_${PROJECT_VERSION}.tar.gz)

add_custom_command(
  OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/test-output.txt
  COMMAND Rscript -e "sink(file='test-output.txt'); devtools::test('waitup'); sink()"
  DEPENDS ${CMAKE_CURRENT_BINARY_DIR}/waitup_${PROJECT_VERSION}.tar.gz
  VERBATIM)

add_custom_target(
  waitup-tests ALL DEPENDS 
    waitup-pkg
    ${CMAKE_CURRENT_BINARY_DIR}/test-output.txt)

