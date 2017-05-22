
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

add_custom_command(
  OUTPUT bugfix-math
  COMMAND git clone https://github.com/stan-dev/cmdstan.git bugfix-cmdstan
  COMMAND git -C bugfix-cmdstan submodule update --init --recursive
  COMMAND git clone https://github.com/stan-dev/math.git bugfix-math
  COMMAND git -C bugfix-math checkout --track origin/bugfix/issue-525-d-gamma-cdf-da
  COMMAND rm -rf bugfix-cmdstan/stan/lib/stan_math
  COMMAND cp -R bugfix-math bugfix-cmdstan/stan/lib/stan_math
  COMMAND make -C bugfix-cmdstan bin/stanc
  VERBATIM)

add_custom_target(bugfix ALL DEPENDS bugfix-math)

add_custom_target(stan-dev ALL DEPENDS bugfix stan-repo stanc-binary cmdstan-repo)

