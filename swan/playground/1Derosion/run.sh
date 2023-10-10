#!/bin/bash
emailid=vaikundaraman@mps.mpg.de
export OMP_NUM_THREADS=48
mail -s 'BEGAN swan 55 playground/1Derosion' $emailid < /dev/null
\time --output=runtime.dat ./2DMC setup.par > run.dat 2>err.dat
cat err.dat runtime.dat > details.dat
#./2DMC setup.par > run.dat 2>err.dat
mail -s 'swan 55 playground/1Derosion done or error' $emailid < details.dat
