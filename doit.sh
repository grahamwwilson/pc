#!/bin/sh

# 6 data jobs
sbatch job_launch.sh DataPA-HPC0A
sbatch job_launch.sh DataPA-HPC0B
sbatch job_launch.sh DataPA-HPC0C
sbatch job_launch.sh DataPA-HPC1
sbatch job_launch.sh DataW-HPC1
sbatch job_launch.sh DataW-HPC2

# 2 MC jobs
sbatch job_launch.sh MCW-HPCA
sbatch job_launch.sh MCW-HPCB

exit
