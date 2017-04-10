#!/bin/bash

echo "functions {" > functions.stan
for f in $1/gamma.stan.part $1/weibull.stan.part \
  $1/weibull-uniform-mixture.stan.part $1/generalized-gamma.stan.part \
  $1/gamma-exp-sum.stan.part $1/gamma-exp-sum-gamma-mix.stan.part
do cat $f >> functions.stan
done
echo "} model {}" >> functions.stan





