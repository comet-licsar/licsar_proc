#!/bin/bash
ddir=$1
if [ -z $ddir ]; then echo "use dir"; exit; fi

for x in `find $ddir -name '*.tif'`; do sz=`ls $x -al | gawk {'print $5'}`; if [ $sz == 0 ]; then echo $x; fi; done
