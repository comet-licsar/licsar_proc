#!/bin/bash

# Plot baselines based on GAMMA bperp file with current interferograms joined
# and colored by % of unwrapped pixels

#################################################################
# Dependencies:
# Requires that a bunch of files have been created before hand
# with the correct naming, etc.
#
# 1) A file generated using GAMMA base_calc: 
#    -> iran_072A_05090_131313_bp_unw.list
# 2) A list file with all of the Scihub acquisition dates
#    -> 072A_05090_131313_bp.list
# 3) A list file with all of the LiCS database acquisition dates
#    -> 072A_05090_131313_scihub.list
 
# Author: Jonathan Weiss
# Date: Sep 2018
#################################################################

usage() { 
  echo " " 1>&2; 
  echo "Purpose: Make a QC baseline plot" 1>&2;
  echo "Usage: baseline_qc_plot.sh [arg1] [arg2] [arg3]" 1>&2; 
  echo " where:" 1>&2;
  echo " [arg1] bperp_file: input Gamma bperp file with unw pixel %" 1>&2;
  echo " [arg2] frame: LiCS frame ID (e.g. 072A_05090_131313) " 1>&2;
  echo " " 1>&2;
  echo " Read header info in script for more info on input files " 1>&2;
  echo " " 1>&2;
  exit 1; 
}

if [ "$#" -ne 2 ]; then
  usage
else
infile=$1
frame=$2
fi

#infile=$1
#location=$2
#frame=$3
title=${frame}
psfile=$title.ps

min_t=`gmtinfo -C $infile | awk '{print $11-50}'`
max_t=`gmtinfo -C $infile | awk '{print $12+100}'`
min_p=`gmtinfo -C $infile | awk '{print $15-150}'`
max_p=`gmtinfo -C $infile | awk '{print $16+50}'`
proj="-JX25/10"
reg="-R2014-09-01T/2019-02-01T/$min_p/$max_p"

#cp /Volumes/a1/insar/sentinel1/COMET_priority_areas/gmt.conf .
#scripts_dir=/group_workspaces/cems2/nceo_geohazards/projects/LiCS/insar_proc/priority_areas/scripts
scripts_dir=/gws/nopw/j04/nceo_geohazards_vol1/projects/LiCS/insar_proc/priority_areas/scripts
cp ${scripts_dir}/gmt.conf `pwd`

gmtset PS_MEDIA a4
gmtset MAP_TICK_LENGTH_SECONDARY 6p/3p
gmtset MAP_LABEL_OFFSET 0.1c

psbasemap $proj $reg -Bpxa3Of1oneg30d -Bsxa1YS -Bpy50+l"Perpendicular Baseline (m)" -BWSen+t"$title" -K -Y5 > $psfile

#Join current IFGs
# IFS=$'\n'
# touch tmp.txt
# for line in $(cat ./$infile)
# do
#   echo "$line" | awk '{print $2,$8}' >> tmp.txt
#   echo "$line" | awk '{print $3,$9}' >> tmp.txt
#   echo ">" >> tmp.txt
#   psxy tmp.txt $proj $reg -O -K -Wgreen >> $psfile
# done

#Join current IFGs and color by % pixels unwrapped
makecpt -Chot -T0/100/.1 -I > unw_pix.cpt
IFS=$'\n'  
for line in `cat $infile`
do
  color=`echo $line | awk '{printf "%.2f\n", $NF}'`
  echo ">"-Z$color > tmp.txt
  echo $line | awk '{print $2,$8}' >> tmp.txt
  echo $line | awk '{print $3,$9}' >> tmp.txt
  psxy tmp.txt $proj $reg -O -K -W.5p -Cunw_pix.cpt  >> $psfile
done


awk '{print $2,$8}' $infile | uniq | psxy $proj $reg -O -K -Sc0.15 -Gblack >> $psfile
awk '{print $3,$9}' $infile | uniq | psxy $proj $reg -O -K -Sc0.15 -Gblack >> $psfile

awk '{print $1,'$min_p'+10}' ${frame}_scihub.list | psxy $proj $reg -O -K -Sc0.15 -Gblue >> $psfile
awk '{print $1,'$min_p'+20}' ${frame}_db.list | psxy $proj $reg -O -K -Sc0.15 -Gred >> $psfile
awk '{print $2,'$min_p'+30}' $infile | uniq | psxy $proj $reg -O -K -Sc0.15 -Gblack >> $psfile
awk '{print $3,'$min_p'+30}' $infile | uniq | psxy $proj $reg -O -K -Sc0.15 -Gblack >> $psfile
#awk '{print 13,'$min_p'+30}' ${location}_${frame}_rslc.list | psxy $proj $reg -O -K -Sc0.15 -Gblack >> $psfile
gmtset MAP_LABEL_OFFSET 0.2c
psscale -Cunw_pix.cpt -D21/2.5/5/0.2h -Ba20f10:"Unwrapped Pixels (%)": -F+gwhite -O >> $psfile
psconvert $psfile -A -P -Tg
