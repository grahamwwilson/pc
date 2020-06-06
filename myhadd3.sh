#!/bin/sh
#
# For new root files
#

date
#. ~/setup.sh

date
# Deal for now with crashing subset PC_data1-3.root by excluding it from the sum
rm PC_data.root
hadd PC_data.root PC_data0.root PC_data1-0.root PC_data1-1.root PC_data1-2.root PC_data1-4.root PC_data1-5.root PC_data1-6.root PC_data1-7.root PC_data1-8.root PC_data1-9.root PC_data2.root

date

rm PC_mc.root
hadd PC_mc.root PC_mc0.root PC_mc1.root PC_mc2.root

date
exit
