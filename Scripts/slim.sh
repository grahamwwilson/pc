#!/bin/sh
#
# Make reduced TTrees using all the .root files in this directory
# (takes 195 minutes for 459 GB running interactively).
#
. ~/setup.sh

ls -1 Run2018_*[!_][!R][!1].root >MyFileList.txt
listfile=MyFileList.txt

date
for file in $(cat ${listfile})

do

   NAME=$(basename ${file} .root)
   echo $NAME
   ~/pc/SlimTree ${NAME}
   date

done

rm MyFileList.txt

exit
