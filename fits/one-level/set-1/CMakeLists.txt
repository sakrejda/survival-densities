
add_custom_command(
  OUTPUT submit.sh
  COMMAND Rscript -e "source('${CMAKE_CURRENT_SOURCE_DIR}/make-fit-targets.R')" ../../../models/stan-lang/full/ ../../../simulation/one-level/set-1/data/ "bsub -q short -W 0:30 -n1 -R 'rusage[mem=512]'" ./submit.sh
  VERBATIM
)

add_custom_target(
  fit-targets ALL DEPENDS submit.sh
)

