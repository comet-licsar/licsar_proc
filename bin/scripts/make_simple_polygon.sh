#!/bin/bash

#################################################################
# hack to make simple (5 coordinate) polygon from LiCS frame coordinates
#
# Dependencies:
# Uses GMT commands so need the following in your path
#
# export PATH=/group_workspaces/cems2/nceo_geohazards/projects/COMET/DavidM/software/local/bin:$PATH
# MANPATH=/group_workspaces/cems2/nceo_geohazards/projects/COMET/jweiss/man:$MANPATH; export MANPATH
#
# Author: Jonathan Weiss
# Date: 2018/07/17
#################################################################

usage() { 
  echo " " 1>&2; 
  echo " Hack to make simple (5 coordinate) polygon from LiCS frame coordinates " 1>&2; 
  echo " Output polygon has a small buffer of 0.01 degrees; " 1>&2; 
  echo "Usage: " 1>&2; 
  echo " Make_simple_polygon.sh LiCS_frame_ID-poly.txt " 1>&2; 
  echo "  " 1>&2; 
  exit 1; 
}

if [ "$#" -ne 1 ]; then
  usage
else
  file=$1
fi

module load gmt

minLon=`gmtinfo $file | sed 's/<//g' | sed 's/>//g' | sed 's/\// /g' | awk '{printf "%5.2f\n",$5-.01}'`
maxLon=`gmtinfo $file | sed 's/<//g' | sed 's/>//g' | sed 's/\// /g' | awk '{printf "%5.2f\n",$6+.01}'`
minLat=`gmtinfo $file | sed 's/<//g' | sed 's/>//g' | sed 's/\// /g' | awk '{printf "%5.2f\n",$7-.01}'`
maxLat=`gmtinfo $file | sed 's/<//g' | sed 's/>//g' | sed 's/\// /g' | awk '{printf "%5.2f\n",$8+.01}'`

outfile=`echo $file | sed 's/-poly.txt/.xy/g'`

echo $minLon $minLat > $outfile
echo $minLon $maxLat >> $outfile
echo $maxLon $maxLat >> $outfile
echo $maxLon $minLat >> $outfile
echo $minLon $minLat >> $outfile
