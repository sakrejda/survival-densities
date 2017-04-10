
file (GLOB target-file-path *.Rnw) 
string (REGEX REPLACE "/.*/" "" target-file ${target-file-path} )
string (REGEX REPLACE "/[^/]+\\.Rnw" "" target-path ${target-file-path} )
string (REGEX REPLACE "\\.Rnw" "" target-name ${target-file} )

add_custom_command(
  OUTPUT ${target-name}.tex
  COMMAND cp ${CMAKE_CURRENT_SOURCE_DIR}/${target-file} ${CMAKE_CURRENT_BINARY_DIR}/
  COMMAND ln -f -n -s ${PROJECT_SOURCE_DIR}/models/stan-lang/full full
  COMMAND ln -f -n -s ${PROJECT_SOURCE_DIR}/models/stan-lang/functions functions
  COMMAND ln -f -n -s ${PROJECT_SOURCE_DIR}/references.bib references.bib
  COMMAND ln -f -n -s ${PROJECT_SOURCE_DIR}/code/waitup/waitup/R R
  COMMAND Rscript -e "knitr::knit('${CMAKE_CURRENT_BINARY_DIR}/${target-name}.Rnw')"
  DEPENDS ${target-file-path}
  VERBATIM
)

add_custom_target(
  ${target-name}-tex ALL
  DEPENDS ${target-name}.tex)

add_custom_command(
  OUTPUT ${target-name}.pdf
  COMMAND pdflatex ${target-name}.tex
  COMMAND biber ${target-name}
  COMMAND pdflatex ${target-name}.tex
  DEPENDS ${target-name}-tex)

add_custom_target(
  ${target-name}-pdf ALL
  DEPENDS ${target-name}.pdf)

