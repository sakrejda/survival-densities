
add_custom_command(
  OUTPUT stannis_1.0.tar.gz
  COMMAND git clone https://github.com/sakrejda/stannis.git
  COMMAND R CMD build stannis
  COMMAND Rscript -e "install.packages(pkgs='stannis_1.0.tar.gz')" 
  VERBATIM)

add_custom_target(
  stannis-build ALL DEPENDS stannis_1.0.tar.gz
)


