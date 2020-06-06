#!/bin/sh

# 3 data jobs
#sbatch job_launch2.sh data0

#sbatch job_launch2.sh data1
sbatch job_launch2.sh data1-30
sbatch job_launch2.sh data1-31
sbatch job_launch2.sh data1-32
sbatch job_launch2.sh data1-33
sbatch job_launch2.sh data1-34
sbatch job_launch2.sh data1-35
sbatch job_launch2.sh data1-36
sbatch job_launch2.sh data1-37
sbatch job_launch2.sh data1-38
sbatch job_launch2.sh data1-39

#sbatch job_launch2.sh data2

# 3 MC jobs
#sbatch job_launch2.sh mc0
#sbatch job_launch2.sh mc1
#sbatch job_launch2.sh mc2

exit
