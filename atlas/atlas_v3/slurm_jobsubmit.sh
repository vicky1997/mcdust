#!/bin/bash
#SBATCH --job-name=2DMC_atlas_v3                   # Job name
#SBATCH --partition=intel                      # Partition to use
#SBATCH --qos=intel_default                    # QOS to be used
#SBATCH --mail-type=BEGIN,END,FAIL,ALL         # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=vaikundaraman@mps.mpg.de   # Where to send mail 
#SBATCH --nodes=1                              # Allocate n nodes
#SBATCH --cpus-per-task=36                     # Number of CPU cores per task
#SBATCH --time=72:00:00                        # Time limit hrs:min:sec
#SBATCH -o run-openMP_%j.out                   # redirect Standard output to file
#SBATCH -e run-openMP_%j.err                   # same for STDERR

# Set OMP_NUM_THREADS to the same value as -c
# with a fallback in case it isn't set.
# SLURM_CPUS_PER_TASK is set to the value of -c, but only if -c is explicitly set
if [ -n "$SLURM_CPUS_PER_TASK" ]; then
  omp_threads=$SLURM_CPUS_PER_TASK
else
  omp_threads=36
fi
export OMP_NUM_THREADS=$omp_threads

./2DMC setup.par
# End of submit file
