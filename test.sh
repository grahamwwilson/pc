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

sbatch job_launch.sh testdata_R${MODE}
sbatch job_launch.sh testmc_R${MODE}

exit
