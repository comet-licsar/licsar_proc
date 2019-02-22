#!/bin/bash
#this script will download a given file (zip file) from alaska facility
URL=https://datapool.asf.alaska.edu/SLC
export WGETRC=$LiCSAR_configpath/alaska_credentials_wget

if [ -z $1 ]; then echo "Parameter is the file name, e.g. S1A....zip"; exit; fi
if [ `echo $1 | grep -c zip` -eq 0 ]; then echo "Parameter is the file name, e.g. S1A....zip"; exit; fi
line=$1

satellite=`echo $line | awk '{print substr($1,1,3)}'`
if [ $satellite == "S1A" ]; then
    satstr='SA'
elif [ $satellite == "S1B" ]; then
    satstr='SB'
fi

wget -c $URL/$satstr/$line
