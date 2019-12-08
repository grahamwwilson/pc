#!/bin/sh
#
# Update to slimmed versions
#

# 6 data jobs
sbatch job_launch.sh DataPA-HPC0A_R1
sbatch job_launch.sh DataPA-HPC0B_R1
sbatch job_launch.sh DataPA-HPC0C_R1
sbatch job_launch.sh DataPAW-HPC1_R1
sbatch job_launch.sh DataW-HPC1_R1
sbatch job_launch.sh DataW-HPC2_R1

# 2 MC jobs
sbatch job_launch.sh MCW-HPCA_R1
sbatch job_launch.sh MCW-HPCB_R1

exit
