#!/bin/sh

listfile=MyFileListNoSuffix.txt

date
for file in $(cat ${listfile})

do

   ~/pc/SlimTree ${file}
   date

done

exit
