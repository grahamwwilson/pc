#!/bin/sh
#
# For new root files
#

date
. ~/setup.sh

date
rm PC_data.root
hadd PC_data.root PC_data0.root PC_data1.root PC_data2.root

date

rm PC_mc.root
hadd PC_mc.root PC_mc0.root PC_mc1.root PC_mc2.root

date
exit
