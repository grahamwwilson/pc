#!/bin/sh

# 3 data jobs
#sbatch job_launch2.sh data0

#sbatch job_launch2.sh data1
sbatch job_launch2.sh data1-390
sbatch job_launch2.sh data1-391
sbatch job_launch2.sh data1-392
sbatch job_launch2.sh data1-393
sbatch job_launch2.sh data1-394
sbatch job_launch2.sh data1-395
sbatch job_launch2.sh data1-396
sbatch job_launch2.sh data1-397
sbatch job_launch2.sh data1-398
sbatch job_launch2.sh data1-399

#sbatch job_launch2.sh data2

# 3 MC jobs
#sbatch job_launch2.sh mc0
#sbatch job_launch2.sh mc1
#sbatch job_launch2.sh mc2

exit
