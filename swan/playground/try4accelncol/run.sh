#!/bin/bash
emailid=vaikundaraman@mps.mpg.de
export OMP_NUM_THREADS=48
mail -s 'BEGAN swan 56 swan/try4/accelncol' $emailid < /dev/null
#\time --output=runtime.out ./2DMC setup.par > run.out 2>err.out
#cat err.out runtime.out > details.out
./2DMC setup.par > run.dat 2>err.dat
mail -s 'swan 56 swan/try4/try4accelncol done or error' $emailid < err.dat
