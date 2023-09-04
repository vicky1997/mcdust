#!/bin/bash
email = vaikundaraman@mps.mpg.de
export OMP_NUM_THREADS=48
mail -s 'run began' $email < /dev/null
time -o runtime.out ./2DMC setup.par > run.out
mail -s 'run either done or error' $email < runtime.out