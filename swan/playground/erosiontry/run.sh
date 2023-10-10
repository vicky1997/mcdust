#!/bin/bash
emailid=vaikundaraman@mps.mpg.de
export OMP_NUM_THREADS=48
mail -s 'BEGAN swan 54 erosion swan/try4/erosiontry' $emailid < /dev/null
#\time --output=runtime.out ./2DMC setup.par > run.out 2>err.out
#cat err.out runtime.out > details.out
./2DMC setup.par > run.dat 2>err.dat
mail -s 'swan 54 swan/try4/erosiontry done or error' $emailid < err.dat
