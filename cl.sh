#!/bin/sh
#
# To compile filename.C with ROOT do
# ./cl.sh filename
# Then the executable can be executed using ./filename
#

#target=$1
target=SlimTree
echo 'Compiling with ROOT libraries '${target}.C

g++ -c PhotonConversionsTree.C `root-config --cflags --glibs`

g++ -g -o ${target} ${target}.C PhotonConversionsTree.o `root-config --cflags --glibs`



exit
