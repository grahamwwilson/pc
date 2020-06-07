#!/bin/sh

INFILE=/home/graham/PC-RootFiles/Data_1396.root

root -l -b -q ${INFILE} -e 'Events->Process("convsel.C")'

exit
