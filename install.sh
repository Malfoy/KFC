#!/bin/bash

folder=$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )
folder+="/bin"
mkdir $folder;
echo "I put binaries in $folder";
( cd thirdparty/smhasher/src; cmake .; make )
make
mv kfc_blue kfc_red kmerCountEvaluator $folder
rm *.o
