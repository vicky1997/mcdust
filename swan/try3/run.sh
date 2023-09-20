#!/bin/bash
emailid=vaikundaraman@mps.mpg.de
touch runtime.out
export OMP_NUM_THREADS=48
mail -s 'BEGAN swan 52 missing particles' $emailid < /dev/null
#\time --output=runtime.out ./2DMC setup.par > run.out 2>err.out
#cat err.out runtime.out > details.out
./2DMC setup.par > run.out 2>err.out
mail -s 'swan 52 missing particles either done or error' $emailid < err.out
