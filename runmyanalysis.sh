#!/bin/sh

date

INPUT=${1-data}

INPUTLIST=${INPUT}.list

#input is: num files, numthreads,  yourdata.list
#0 is all files, optimal threads is 8
#python runmacro.py 0 8 "bigmcdata.list"
python2 runmacro.py 0 8 ${INPUTLIST}

mv Outfile.root PC_${INPUT}.root

date

exit
