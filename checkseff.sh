#!/bin/sh
#
# Get batch job statistics for a sequence of completed jobs
#

jstart=$1
jend=$2

echo "Checking batch job statistics for job "${jstart}" to "${jend}
#echo "to "${jend}

i=${jstart}

while [ $i -le ${jend} ]
do
#    echo "job " $i
    seff $i
    echo " " 
    let "i+=1"
done

exit
