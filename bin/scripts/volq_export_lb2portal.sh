#!/bin/bash

if [ -z $1 ]; then
 echo "Usage: run in folder with TS directory containing cum*h5"
 echo "(will run with any one param, e.g. 1)"
 exit
fi

hdf=TS*/cum_filt.h5
LiCSBAS_out2nc.py -i $hdf -o temp.nc --cf
nccopy -4 -d 4 -s -u temp.nc out.64.nc
