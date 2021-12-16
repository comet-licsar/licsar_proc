#!/bin/bash

source $LiCSARpath/lib/LiCSAR_bash_lib.sh


if [ -z $1 ]; then
  cat << End_of_Usage

  Purpose: Fix geocoding error in some frame geotiffs (will shift it w.r.t. hgt)


End_of_Usage
  exit
fi

frame=$1
track=`track_from_frame $frame`

cdpub $frame
hgt=`pwd`/metadata/$frame.geo.hgt.tif

coordhgt=`gdalinfo $hgt | grep Center | gawk {'print $3'} | cut -d ',' -f1`

for tif in `ls metadata/*tif`; do
 coord=`gdalinfo $tif | grep Center | gawk {'print $3'} | cut -d ',' -f1`
 if [ ! $coordhgt == $coord ]; then
   echo Correcting $tif
   gdalwarp2match.py $tif $hgt temp_warped.tif >/dev/null 2>/dev/null
   rm -f $tif
   gdal_translate -of GTiff -co COMPRESS=DEFLATE temp_warped.tif $tif >/dev/null 2>/dev/null
   rm temp_warped.tif
 fi
done

for epoch in `ls epochs | grep 20`; do
 for tif in `ls epochs/$epoch/*tif`; do
  coord=`gdalinfo $tif | grep Center | gawk {'print $3'} | cut -d ',' -f1`
  if [ ! $coordhgt == $coord ]; then
    echo Correcting $tif
    gdalwarp2match.py $tif $hgt temp_warped.tif >/dev/null 2>/dev/null
    rm -f $tif
    gdal_translate -of GTiff -co COMPRESS=DEFLATE temp_warped.tif $tif >/dev/null 2>/dev/null
    rm temp_warped.tif
  fi
 done
done


for pair in `ls interferograms | grep 20`; do
 for tif in `ls interferograms/$pair/*tif`; do
   coord=`gdalinfo $tif | grep Center | gawk {'print $3'} | cut -d ',' -f1`
   if [ ! $coordhgt == $coord ]; then
     echo Correcting $tif
     gdalwarp2match.py $tif $hgt temp_warped.tif >/dev/null 2>/dev/null
     rm -f $tif
     gdal_translate -of GTiff -co COMPRESS=DEFLATE temp_warped.tif $tif >/dev/null 2>/dev/null
     rm temp_warped.tif
   else
     echo "this is ok: "$tif
   fi
 done
done
