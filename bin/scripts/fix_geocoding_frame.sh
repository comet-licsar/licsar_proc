#!/bin/bash

source $LiCSARpath/lib/LiCSAR_bash_lib.sh


if [ -z $1 ]; then
  cat << End_of_Usage

  Purpose: Fix geocoding error in some frame geotiffs (will shift it w.r.t. hgt)


End_of_Usage
  exit
fi

# to make smaller:
# fr1=143A_04800_131313
# cdpub $fr1;cd metadata
# gdalwarp -tr 0.001 0.001 $fr1.geo.hgt.tif $fr1.geo.hgt.tif.2.tif
# rm $fr1.geo.hgt.tif
# gdal_translate -of GTiff -co COMPRESS=DEFLATE $fr1.geo.hgt.tif.2.tif $fr1.geo.hgt.tif
# rm $fr1.geo.hgt.tif.2.tif
frame=$1
track=`track_from_frame $frame`

cdpub $frame
hgt=`pwd`/metadata/$frame.geo.hgt.tif

coordhgt=`gdalinfo $hgt | grep Center | gawk {'print $3'} | cut -d ',' -f1`
sizehgt=`gdalinfo $hgt | grep Size | head -n 1 | sed 's/ //g'`

for tif in `ls metadata/*tif`; do
 coord=`gdalinfo $tif | grep Center | gawk {'print $3'} | cut -d ',' -f1`
 if [ ! $coordhgt == $coord ]; then
   corit=1
 else
   sizetif=`gdalinfo $tif | grep Size | head -n 1 | sed 's/ //g'`
   if [ ! $sizehgt == $sizetif ]; then
     corit=1
   else
     corit=0
   fi
 fi
 if [ $corit == 1 ]; then
   echo Correcting $tif
   gdalwarp2match.py $tif $hgt temp_warped.tif >/dev/null 2>/dev/null
   rm -f $tif
   gdal_translate -of GTiff -co COMPRESS=DEFLATE temp_warped.tif $tif >/dev/null 2>/dev/null
   rm temp_warped.tif
 else
  echo "this is ok: "$tif
 fi
done

for epoch in `ls epochs | grep 20`; do
 for tif in `ls epochs/$epoch/*tif`; do
 coord=`gdalinfo $tif | grep Center | gawk {'print $3'} | cut -d ',' -f1`
 if [ ! $coordhgt == $coord ]; then
   corit=1
 else
   sizetif=`gdalinfo $tif | grep Size | head -n 1 | sed 's/ //g'`
   if [ ! $sizehgt == $sizetif ]; then
     corit=1
   else
     corit=0
   fi
 fi
 if [ $corit == 1 ]; then
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


for pair in `ls interferograms | grep 20`; do
 for tif in `ls interferograms/$pair/*tif`; do
 coord=`gdalinfo $tif | grep Center | gawk {'print $3'} | cut -d ',' -f1`
 if [ ! $coordhgt == $coord ]; then
   corit=1
 else
   sizetif=`gdalinfo $tif | grep Size | head -n 1 | sed 's/ //g'`
   if [ ! $sizehgt == $sizetif ]; then
     corit=1
   else
     corit=0
   fi
 fi
 if [ $corit == 1 ]; then
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
