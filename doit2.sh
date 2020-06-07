#!/bin/sh

# 3 data jobs
#sbatch job_launch2.sh data0

sbatch job_launch2.sh data1
#sbatch job_launch2.sh data1-0
#sbatch job_launch2.sh data1-1
#sbatch job_launch2.sh data1-2
#sbatch job_launch2.sh data1-3
#sbatch job_launch2.sh data1-4
#sbatch job_launch2.sh data1-5
#sbatch job_launch2.sh data1-6
#sbatch job_launch2.sh data1-7
#sbatch job_launch2.sh data1-8
#sbatch job_launch2.sh data1-9

#sbatch job_launch2.sh data2

# 3 MC jobs
#sbatch job_launch2.sh mc0
#sbatch job_launch2.sh mc1
#sbatch job_launch2.sh mc2

exit
