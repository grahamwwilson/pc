#!/bin/sh

. ~/setup.sh

date
hadd PC_MCHPC.root PC_MCW-HPCA.root PC_MCW-HPCB.root
date
hadd PC_DataHPC.root \
     PC_DataPA-HPC0A.root PC_DataPA-HPC0B.root PC_DataPA-HPC0C.root \
     PC_DataPA-HPC1.root PC_DataW-HPC1.root PC_DataW-HPC2.root
date
exit
