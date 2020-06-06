#!/bin/sh

# 3 data jobs
sbatch job_launch2.sh data0
sbatch job_launch2.sh data1
sbatch job_launch2.sh data2

# 3 MC jobs
sbatch job_launch2.sh mc0
sbatch job_launch2.sh mc1
sbatch job_launch2.sh mc2

exit
