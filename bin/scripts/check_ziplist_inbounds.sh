#!/bin/bash

infile=$1
polygon=$2

for zipfile1 in `cat $infile`; do
  inbounds=`zipSAFE_inbounds_frame.py $zipfile1 $polygon`
  #echo $inbounds $i
  if [ "$inbounds" == "1" ]; then
    echo $zipfile1
  fi
done
