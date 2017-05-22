
add_custom_command(
  OUTPUT stannis_1.0.tar.gz
  COMMAND rm -rf ${CMAKE_CURRENT_BINARY_DIR}/stannis
  COMMAND git clone https://github.com/sakrejda/stannis.git
  COMMAND R_LIBS_USER=${R_LIB_DIR} R CMD build stannis
  COMMAND Rscript -e "install.packages(pkgs='stannis_1.0.tar.gz')" 
  COMMAND R CMD INSTALL ${CMAKE_CURRENT_BINARY_DIR}/stannis_1.0.tar.gz --library=${R_LIB_DIR} 
  VERBATIM)

add_custom_target(
  stannis-build ALL DEPENDS stannis_1.0.tar.gz
)


