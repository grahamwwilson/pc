#!/bin/sh
#
# By default ths now uses slimmed versions (MODE=1)      R1
#
# Specify MODE=0 to get the unslimmed version            R0
# (this may still need some changes to the code - to check) 
#
MODE=${1-1}
echo 'MODE set to '
echo ${MODE}

# 7 data jobs
sbatch job_launch.sh DataPA-HPC0A_R${MODE}
sbatch job_launch.sh DataPA-HPC0B_R${MODE}
sbatch job_launch.sh DataPA-HPC0C_R${MODE}
sbatch job_launch.sh DataPAW-HPC1A_R${MODE}
sbatch job_launch.sh DataPAW-HPC1B_R${MODE}
sbatch job_launch.sh DataW-HPC1_R${MODE}
sbatch job_launch.sh DataW-HPC2_R${MODE}

# 2 MC jobs
sbatch job_launch.sh MCW-HPCA_R${MODE}
sbatch job_launch.sh MCW-HPCB_R${MODE}

exit
