
add_custom_command(
  OUTPUT gamma-p1.stan gamma-p1
  COMMAND cat ${CMAKE_CURRENT_SOURCE_DIR}/head.stan.part > gamma-p1.stan
  COMMAND cat ${CMAKE_CURRENT_SOURCE_DIR}/../functions/gamma.stan.part >> gamma-p1.stan
  COMMAND cat ${CMAKE_CURRENT_SOURCE_DIR}/tail.stan.part >> gamma-p1.stan
  COMMAND cat ${CMAKE_CURRENT_SOURCE_DIR}/../data/truncated-delay.stan.part >> gamma-p1.stan
  COMMAND cat ${CMAKE_CURRENT_SOURCE_DIR}/../bodies/gamma-p1.stan.part >> gamma-p1.stan
  COMMAND make -C ${CMDSTAN_FIX_DIR} ${CMAKE_CURRENT_BINARY_DIR}/gamma-p1
  DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/../functions/gamma.stan.part 
          ${CMAKE_CURRENT_SOURCE_DIR}/../data/truncated-delay.stan.part 
          ${CMAKE_CURRENT_SOURCE_DIR}/../bodies/gamma-p1.stan.part
  VERBATIM)

add_custom_command(
  OUTPUT weibull-p1.stan weibull-p1
  COMMAND cat ${CMAKE_CURRENT_SOURCE_DIR}/head.stan.part > weibull-p1.stan
  COMMAND cat ${CMAKE_CURRENT_SOURCE_DIR}/../functions/weibull.stan.part >> weibull-p1.stan
  COMMAND cat ${CMAKE_CURRENT_SOURCE_DIR}/tail.stan.part >> weibull-p1.stan
  COMMAND cat ${CMAKE_CURRENT_SOURCE_DIR}/../data/truncated-delay.stan.part >> weibull-p1.stan
  COMMAND cat ${CMAKE_CURRENT_SOURCE_DIR}/../bodies/weibull-p1.stan.part >> weibull-p1.stan
  COMMAND make -C ${CMDSTAN_FIX_DIR} ${CMAKE_CURRENT_BINARY_DIR}/weibull-p1
  DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/../functions/weibull.stan.part 
          ${CMAKE_CURRENT_SOURCE_DIR}/../data/truncated-delay.stan.part 
          ${CMAKE_CURRENT_SOURCE_DIR}/../bodies/weibull-p1.stan.part
  VERBATIM)

add_custom_command(
  OUTPUT weibull-uniform-mixture-p1.stan weibull-uniform-mixture-p1
  COMMAND cat ${CMAKE_CURRENT_SOURCE_DIR}/head.stan.part > weibull-uniform-mixture-p1.stan
  COMMAND cat ${CMAKE_CURRENT_SOURCE_DIR}/../functions/weibull.stan.part >> weibull-uniform-mixture-p1.stan
  COMMAND cat ${CMAKE_CURRENT_SOURCE_DIR}/../functions/weibull-uniform-mixture.stan.part >> weibull-uniform-mixture-p1.stan
  COMMAND cat ${CMAKE_CURRENT_SOURCE_DIR}/tail.stan.part >> weibull-uniform-mixture-p1.stan
  COMMAND cat ${CMAKE_CURRENT_SOURCE_DIR}/../data/truncated-delay.stan.part >> weibull-uniform-mixture-p1.stan
  COMMAND cat ${CMAKE_CURRENT_SOURCE_DIR}/../bodies/weibull-uniform-mixture-p1.stan.part >> weibull-uniform-mixture-p1.stan
  COMMAND make -C ${CMDSTAN_FIX_DIR} ${CMAKE_CURRENT_BINARY_DIR}/weibull-uniform-mixture-p1
  DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/../functions/weibull-uniform-mixture.stan.part 
          ${CMAKE_CURRENT_SOURCE_DIR}/../data/truncated-delay.stan.part 
          ${CMAKE_CURRENT_SOURCE_DIR}/../bodies/weibull-uniform-mixture-p1.stan.part
  VERBATIM)

add_custom_command(
  OUTPUT generalized-gamma-p1.stan generalized-gamma-p1
  COMMAND cat ${CMAKE_CURRENT_SOURCE_DIR}/head.stan.part > generalized-gamma-p1.stan
  COMMAND cat ${CMAKE_CURRENT_SOURCE_DIR}/../functions/generalized-gamma.stan.part >> generalized-gamma-p1.stan
  COMMAND cat ${CMAKE_CURRENT_SOURCE_DIR}/tail.stan.part >> generalized-gamma-p1.stan
  COMMAND cat ${CMAKE_CURRENT_SOURCE_DIR}/../data/truncated-delay.stan.part >> generalized-gamma-p1.stan
  COMMAND cat ${CMAKE_CURRENT_SOURCE_DIR}/../bodies/generalized-gamma-p1.stan.part >> generalized-gamma-p1.stan
  COMMAND make -C ${CMDSTAN_FIX_DIR} ${CMAKE_CURRENT_BINARY_DIR}/generalized-gamma-p1
  DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/../functions/generalized-gamma.stan.part 
          ${CMAKE_CURRENT_SOURCE_DIR}/../data/truncated-delay.stan.part 
          ${CMAKE_CURRENT_SOURCE_DIR}/../bodies/generalized-gamma-p1.stan.part
  VERBATIM)

add_custom_command(
  OUTPUT generalized-gamma-p2.stan generalized-gamma-p2
  COMMAND cat ${CMAKE_CURRENT_SOURCE_DIR}/head.stan.part > generalized-gamma-p2.stan
  COMMAND cat ${CMAKE_CURRENT_SOURCE_DIR}/../functions/generalized-gamma.stan.part >> generalized-gamma-p2.stan
  COMMAND cat ${CMAKE_CURRENT_SOURCE_DIR}/tail.stan.part >> generalized-gamma-p2.stan
  COMMAND cat ${CMAKE_CURRENT_SOURCE_DIR}/../data/truncated-delay.stan.part >> generalized-gamma-p2.stan
  COMMAND cat ${CMAKE_CURRENT_SOURCE_DIR}/../bodies/generalized-gamma-p2.stan.part >> generalized-gamma-p2.stan
  COMMAND make -C ${CMDSTAN_FIX_DIR} ${CMAKE_CURRENT_BINARY_DIR}/generalized-gamma-p2
  DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/../functions/generalized-gamma.stan.part 
          ${CMAKE_CURRENT_SOURCE_DIR}/../data/truncated-delay.stan.part 
          ${CMAKE_CURRENT_SOURCE_DIR}/../bodies/generalized-gamma-p2.stan.part
  VERBATIM)

add_custom_command(
  OUTPUT gamma-exp-sum-p1.stan gamma-exp-sum-p1
  COMMAND cat ${CMAKE_CURRENT_SOURCE_DIR}/head.stan.part > gamma-exp-sum-p1.stan
  COMMAND cat ${CMAKE_CURRENT_SOURCE_DIR}/../functions/gamma-exp-sum.stan.part >> gamma-exp-sum-p1.stan
  COMMAND cat ${CMAKE_CURRENT_SOURCE_DIR}/tail.stan.part >> gamma-exp-sum-p1.stan
  COMMAND cat ${CMAKE_CURRENT_SOURCE_DIR}/../data/truncated-delay.stan.part >> gamma-exp-sum-p1.stan
  COMMAND cat ${CMAKE_CURRENT_SOURCE_DIR}/../bodies/gamma-exp-sum-p1.stan.part >> gamma-exp-sum-p1.stan
  COMMAND make -C ${CMDSTAN_FIX_DIR} ${CMAKE_CURRENT_BINARY_DIR}/gamma-exp-sum-p1
  DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/../functions/gamma-exp-sum.stan.part 
          ${CMAKE_CURRENT_SOURCE_DIR}/../data/truncated-delay.stan.part 
          ${CMAKE_CURRENT_SOURCE_DIR}/../bodies/gamma-exp-sum-p1.stan.part
  VERBATIM)

add_custom_command(
  OUTPUT gamma-exp-sum-gamma-mix-p1.stan gamma-exp-sum-gamma-mix-p1
  COMMAND cat ${CMAKE_CURRENT_SOURCE_DIR}/head.stan.part > gamma-exp-sum-gamma-mix-p1.stan
  COMMAND cat ${CMAKE_CURRENT_SOURCE_DIR}/../functions/gamma.stan.part >> gamma-exp-sum-gamma-mix-p1.stan
  COMMAND cat ${CMAKE_CURRENT_SOURCE_DIR}/../functions/gamma-exp-sum.stan.part >> gamma-exp-sum-gamma-mix-p1.stan
  COMMAND cat ${CMAKE_CURRENT_SOURCE_DIR}/../functions/gamma-exp-sum-gamma-mix.stan.part >> gamma-exp-sum-gamma-mix-p1.stan
  COMMAND cat ${CMAKE_CURRENT_SOURCE_DIR}/tail.stan.part >> gamma-exp-sum-gamma-mix-p1.stan
  COMMAND cat ${CMAKE_CURRENT_SOURCE_DIR}/../data/truncated-delay.stan.part >> gamma-exp-sum-gamma-mix-p1.stan
  COMMAND cat ${CMAKE_CURRENT_SOURCE_DIR}/../bodies/gamma-exp-sum-gamma-mix-p1.stan.part >> gamma-exp-sum-gamma-mix-p1.stan
  COMMAND make -C ${CMDSTAN_FIX_DIR} ${CMAKE_CURRENT_BINARY_DIR}/gamma-exp-sum-gamma-mix-p1
  DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/../functions/gamma-exp-sum-gamma-mix.stan.part 
          ${CMAKE_CURRENT_SOURCE_DIR}/../data/truncated-delay.stan.part 
          ${CMAKE_CURRENT_SOURCE_DIR}/../bodies/gamma-exp-sum-gamma-mix-p1.stan.part
  VERBATIM)

add_custom_command(
  OUTPUT gamma-exp-sum-gamma-mix-p2.stan gamma-exp-sum-gamma-mix-p2
  COMMAND cat ${CMAKE_CURRENT_SOURCE_DIR}/head.stan.part > gamma-exp-sum-gamma-mix-p2.stan
  COMMAND cat ${CMAKE_CURRENT_SOURCE_DIR}/../functions/gamma.stan.part >> gamma-exp-sum-gamma-mix-p2.stan
  COMMAND cat ${CMAKE_CURRENT_SOURCE_DIR}/../functions/gamma-exp-sum.stan.part >> gamma-exp-sum-gamma-mix-p2.stan
  COMMAND cat ${CMAKE_CURRENT_SOURCE_DIR}/../functions/gamma-exp-sum-gamma-mix.stan.part >> gamma-exp-sum-gamma-mix-p2.stan
  COMMAND cat ${CMAKE_CURRENT_SOURCE_DIR}/tail.stan.part >> gamma-exp-sum-gamma-mix-p2.stan
  COMMAND cat ${CMAKE_CURRENT_SOURCE_DIR}/../data/truncated-delay.stan.part >> gamma-exp-sum-gamma-mix-p2.stan
  COMMAND cat ${CMAKE_CURRENT_SOURCE_DIR}/../bodies/gamma-exp-sum-gamma-mix-p2.stan.part >> gamma-exp-sum-gamma-mix-p2.stan
  COMMAND make -C ${CMDSTAN_FIX_DIR} ${CMAKE_CURRENT_BINARY_DIR}/gamma-exp-sum-gamma-mix-p2
  DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/../functions/gamma-exp-sum-gamma-mix.stan.part 
          ${CMAKE_CURRENT_SOURCE_DIR}/../functions/gamma-exp-sum.stan.part 
          ${CMAKE_CURRENT_SOURCE_DIR}/../functions/gamma.stan.part 
          ${CMAKE_CURRENT_SOURCE_DIR}/../data/truncated-delay.stan.part 
          ${CMAKE_CURRENT_SOURCE_DIR}/../bodies/gamma-exp-sum-gamma-mix-p2.stan.part
  VERBATIM)

add_custom_target(compiled-gamma-p1-fix ALL DEPENDS 
  gamma-p1.stan gamma-p1)
add_custom_target(compiled-weibull-p1-fix ALL DEPENDS 
  weibull-p1.stan weibull-p1)
add_custom_target(compiled-weibull-uniform-mixture-p1-fix ALL DEPENDS 
  weibull-uniform-mixture-p1.stan weibull-uniform-mixture-p1)
add_custom_target(compiled-generalized-gamma-p1-fix ALL DEPENDS 
  generalized-gamma-p1.stan generalized-gamma-p1)
add_custom_target(compiled-generalized-gamma-p2-fix ALL DEPENDS 
  generalized-gamma-p2.stan generalized-gamma-p2)
add_custom_target(compiled-gamma-exp-sum-p1-fix ALL DEPENDS 
  gamma-exp-sum-p1.stan gamma-exp-sum-p1)
add_custom_target(compiled-gamma-exp-sum-gamma-mix-p1-fix ALL DEPENDS 
  gamma-exp-sum-gamma-mix-p1.stan gamma-exp-sum-gamma-mix-p1)
add_custom_target(compiled-gamma-exp-sum-gamma-mix-p2-fix ALL DEPENDS 
  gamma-exp-sum-gamma-mix-p2.stan gamma-exp-sum-gamma-mix-p2)


