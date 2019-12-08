#!/bin/sh
# Update with slimmed versions

date
. ~/setup.sh

date
rm PC_MCHPC.root
hadd PC_MCHPC.root PC_MCW-HPCA_R1.root PC_MCW-HPCB_R1.root

date
rm PC_DataHPC.root
hadd PC_DataHPC.root \
     PC_DataPA-HPC0A_R1.root PC_DataPA-HPC0B_R1.root PC_DataPA-HPC0C_R1.root \
     PC_DataPAW-HPC1_R1.root PC_DataW-HPC1_R1.root PC_DataW-HPC2_R1.root

date
exit
