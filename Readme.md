To compile the code: make

parameter file: setup.par

To run the code: 

export OMP_NUM_THREADS=number of CPUs

./2DMC setup.par

or use the script, (you will have to modify the script a bit to add your emails) 

./run.sh

To run the code via a SLURM batch script 

sbatch slurm_jobsubmit.sh

(You need to modify the contents of the script to suit your slurm account and cluster)
