#!/bin/bash

ns=('2' '10' '100')
priors=('flat' 'weak')
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
  	sbatch triton/misspecified_run.sh $n $prior $iter
    done
  done
done
