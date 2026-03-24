#!/bin/bash
#SBATCH -c 1                               # Request one core
#SBATCH -t 1-13:00                         # Runtime in D-HH:MM format
#SBATCH -p medium                          # Partition to run in
#SBATCH --mem-per-cpu=49G                 # Memory per core
#SBATCH -o ElasticNet_%A_%a.out                 # File to which STDOUT will be written, including job ID (%j)
#SBATCH -e ElasticNet_%A_%a.err                 # File to which STDERR will be written, including job ID (%j) 
#SBATCH --array=1-20          # Run array for indexes 1-20, effectively running the command 20x 
module load gcc/9.2.0 R/4.1.2

unset SLURM_CPU_BIND
srun -n 1 Rscript /home/gb178/tuberculosis/MultipleLogisticRegressionWithLasso.R
sleep 5                                 # wait for slurm to get the job status into its database
sacct --format=JobID,Submit,Start,End,State,Partition,ReqTRES%30,CPUTime,MaxRSS,NodeList%30 --units=M -j $SLURM_JOBID
