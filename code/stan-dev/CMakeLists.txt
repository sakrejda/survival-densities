
add_custom_command(
  OUTPUT cmdstan
  COMMAND git clone https://github.com/stan-dev/cmdstan.git
  COMMAND git -C cmdstan submodule update --init --recursive
  VERBATIM)

add_custom_target(
  cmdstan-repo ALL DEPENDS cmdstan)

add_custom_command(
  OUTPUT cmdstan/bin/stanc
  COMMAND make -C cmdstan bin/stanc
  VERBATIM)

add_custom_target(
  stanc-binary ALL DEPENDS cmdstan/bin/stanc cmdstan-repo)

add_custom_command(
  OUTPUT stan
  COMMAND git clone https://github.com/stan-dev/stan.git
  COMMAND git -C stan submodule update --init --recursive
  VERBATIM)

add_custom_target(
  stan-repo ALL DEPENDS stan)



