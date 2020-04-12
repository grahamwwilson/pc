#!/bin/sh
# Update with slimmed versions

#
# By default ths now uses slimmed versions (MODE=1)      R1
#
# Specify MODE=0 to get the unslimmed version            R0
# (this may still need some changes to the code - to check) 
#
MODE=${1-1}
echo 'MODE set to '
echo ${MODE}

date
. ~/setup.sh

date
rm PC_MCHPC_R${MODE}.root
hadd PC_MCHPC_R${MODE}.root PC_MCW-HPCA_R${MODE}.root PC_MCW-HPCB_R${MODE}.root

date
rm PC_DataHPC_R${MODE}.root
hadd PC_DataHPC_R${MODE}.root \
     PC_DataPA-HPC0A_R${MODE}.root PC_DataPA-HPC0B_R${MODE}.root PC_DataPA-HPC0C_R${MODE}.root \
     PC_DataPAW-HPC1A_R${MODE}.root PC_DataPAW-HPC1B_R${MODE}.root \
     PC_DataW-HPC1_R${MODE}.root PC_DataW-HPC2_R${MODE}.root

date
exit
