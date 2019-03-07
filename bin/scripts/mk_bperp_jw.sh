#!/bin/bash
#################################################################
# Script to process TOPS S1 images using GAMMA software
# Back to back from SLC to geocoded unwrapped interferograms
#
# Author: Pablo J. Gonzalez (p.j.gonzalez@leeds.ac.uk)
# Date: 2015/10/26
#
# Modified by Jonathan Weiss to  create baseline pair lists that can be
# be parsed and fed to python LiCSAR_03_mk_ifgs.py so if you want to 
# make longer temporal baseline pairs, etc.
#
# Note that script creates file of "good" RSLCs based on the existence
# of the existing of the yyyymmdd.rslc.par file 
#
# Adjust bplim and btlim and examine output bperp_file 
# 
# Date July 2018
#################################################################

usage() { 
  echo " " 1>&2; 
  echo "Purpose: Make a baseline plot" 1>&2;
  echo " " 1>&2;
  exit 1; 
}

if [ "$#" -ne 2 ]; then
  usage
else
  bplim=$1
  btlim=$2
fi

master=`ls geo/????????.hgt | cut -c 5-12`         # Master date
#slcs_list=$1      # List of images to be processed 
ls RSLC/*/????????.rslc.par | sed 's/\// /g' | awk '{print $3}' | cut -c 1-8 | awk '{print $1,"1 0"}' > rslc.list
rlks=20           # Multilooking factor for range
azlks=4          # Multilooking factor for azimuth
outres=0.001         # Output resolution of final products (decimal degrees)
# In case we cropped based on geographical coordinates we have to ensure that we always process the same subswaths 
# That happens when the SLC_tab files are created
if [ -e minmaxIW ]; then 
  IW1=`awk '{print $1}' minmaxIW`
  IW3=`awk '{print $2}' minmaxIW`
else
  IW1=1 
  IW3=3
fi
dem=DEM/dem_crop.dem
dempar=${dem}_par
# Select Bperp (m) and Bt (days) limits
#bplim=$1   # Basically accept all perpendicular baselines   
#btlim=$2    # 365 days temporal baselines
delta_n_max=20 # pairs with max distance of 3 SLC regardless of temporal baseline

#################################################################
#################################################################
#################################################################

#----------------------------------------------------------------------------------------#
# Extract information about the master file (dimension, coordinates, etc.)
if [ -e SLC/${master}/${master}.slc.mli.par ]; then
  width=`awk '$1 == "range_samples:" {print $2}' SLC/${master}/${master}.slc.mli.par`; #echo $width
  length=`awk '$1 == "azimuth_lines:" {print $2}' SLC/${master}/${master}.slc.mli.par`; #echo $length
  reducfac=`echo $width | awk '{if(int($1/1000) > 1) print int($1/1000); else print 1}'`
fi

echo " Running doBperp_Compute step "
echo "   check doBperp_Compute.log if something goes wrong "
logfile=00_doBperp_Compute.log
rm -f $logfile
rm -f SLC_mos*
#for image in `cat $slcs_list | awk '{print $1}'`; do #for image in `ls -d RSLC/2* | sed 's/.*\///'`; do
while read line; do
  image=`echo $line | awk '{print $1}'`; 
  validity=`echo $line | awk '{print $2}'`; 
  if [ "$validity" == "1" ]; then
    echo "RSLC/${image}/${image}.rslc RSLC/${image}/${image}.rslc.par" >> SLC_mosaics_tab_temp
  fi
done < rslc.list 
ls -lhH `awk '{print $2}' SLC_mosaics_tab_temp` > trash; 
grep K trash | awk '{print $NF}' | cut -c 6-13 > trash1
for rslc in `cat trash1`; do grep $rslc SLC_mosaics_tab_temp >> SLC_mosaics_tab; done
\rm trash*

base_calc SLC_mosaics_tab RSLC/$master/$master.rslc.par bperp_file itab 1 0 0 $bplim 1 $btlim $delta_n_max >> $logfile
awk '{print $2,$3}' bperp_file > bperp_ifgs.list
rm -f base.out base_calc.log SLC_mosaics_tab itab
