#!/bin/sh

listfile=MyFileList.txt

for file in $(cat ${listfile})

do

   basename ${file} .root

done

exit
