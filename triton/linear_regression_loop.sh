#!/bin/bash

#ns=('50' '100' '500')
#priors=('flat' 'weak')
ns=('50')
priors=('weak')
n_iters=1000

# loop over n
for n in "${ns[@]}"
do
  # loop over priors
  for prior in "${priors[@]}"
  do
    # loop over iterations
    for iter in $(seq 1 $n_iters)
    do
  	sbatch triton/linear_regression_run.sh $n $prior $iter
    done
  done
done
