#!/bin/sh
#
# Make time-stamped backup of the relevant ROOT files.
# 
DATE=`date "+%b-%d-%Y_%T-%Z"`
echo $DATE
BDIR=Backup_${DATE}
mkdir $BDIR
cp -p PC_data.root $BDIR
cp -p PC_mc.root $BDIR
exit
