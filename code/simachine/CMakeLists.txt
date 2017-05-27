
add_custom_command(
  OUTPUT simachine_1.0.tar.gz
  COMMAND rm -rf simachine
  COMMAND git clone https://github.com/sakrejda/simachine.git
  COMMAND R_LIBS_USER=${R_LIB_DIR} R CMD build simachine
  COMMAND mkdir -p ${R_LIB_DIR}
  COMMAND R_LIBS_USER=${R_LIB_DIR} Rscript -e "install.packages('digest')"
  COMMAND R CMD INSTALL simachine_1.0.tar.gz --library=${R_LIB_DIR}
  DEPENDS r-environment
  VERBATIM)

add_custom_target(
  simachine-build ALL DEPENDS simachine_1.0.tar.gz 
)


