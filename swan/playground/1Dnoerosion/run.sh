#!/bin/bash
emailid=vaikundaraman@mps.mpg.de
export OMP_NUM_THREADS=48
mail -s 'BEGAN swan 54 playground/1Dnoerosion' $emailid < /dev/null
\time --output=runtime.dat ./2DMC setup.par > run.dat 2>err.dat
cat err.dat runtime.dat > details.dat
#./2DMC setup.par > run.dat 2>err.dat
mail -s 'swan 54 playground/1Dnoerosion done or error' $emailid < details.dat
