#!/bin/bash
#this script will download a given file (zip file) from ASF to LiCSAR_SLC and ingest it to db

if [ -z $1 ]; then echo "Parameter is the file name, e.g. S1A....zip"; exit; fi
if [ `echo $1 | grep -c zip` -eq 0 ]; then echo "Parameter is the file name, e.g. S1A....zip"; exit; fi
zipf=$1

cd $LiCSAR_SLC
wget_alaska $zipf
arch2DB.py -f `pwd`/$zipf
cd -
