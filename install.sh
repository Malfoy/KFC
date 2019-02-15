#!/bin/bash

folder=$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )
folder+="/bin"
mkdir $folder;
echo "I put binaries in $folder";
( cd smhasher/src; cmake .; make )
make
mv kfc $folder
rm *.o
