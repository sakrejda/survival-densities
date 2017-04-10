
add_custom_command(
  OUTPUT simachine_1.0.tar.gz
  COMMAND git clone https://github.com/sakrejda/simachine.git
  COMMAND R CMD build simachine
  COMMAND Rscript -e "install.packages(pkgs='simachine_1.0.tar.gz')" 
  VERBATIM)

add_custom_target(
  simachine-build ALL DEPENDS simachine_1.0.tar.gz
)


