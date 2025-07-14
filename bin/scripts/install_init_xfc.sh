#!/bin/bash
if [ `hostname | grep -c vm` -lt 1 ]; then
  echo "please use sci-vm for this script"
  exit
fi
echo "initialising xfc disk"
xfc init
xfcpath=`xfc path`
if [ -z $xfcpath ]; then echo "some error - try running xfc init in terminal"; exit; fi
echo "successfully created in "$xfcpath
echo "export XFCPATH="$xfcpath >> ~/.bashrc
echo "setting custom autodownload dir"
mkdir $xfcpath/SLC
echo $xfcpath/SLC >> $LiCSAR_configpath/autodownloaddirs