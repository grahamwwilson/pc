#!/bin/sh

date

INPUTLIST=${1-data.list}

#input is: num files, numthreads,  yourdata.list
#0 is all files, optimal threads is 8
#python runmacro.py 0 8 "bigmcdata.list"
python runmacro.py 0 8 ${INPUTLIST}

date

exit
