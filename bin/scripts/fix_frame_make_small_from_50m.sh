#!/bin/bash

if [ -z $1 ]; then
  cat << End_of_Usage

  This will process ALL geotiffs decreasing size to 100 m


End_of_Usage
  exit
fi

frame=$1
track=`track_from_frame $frame`

cdpub $frame
hgt=`pwd`/metadata/$frame.geo.hgt.tif
#sz=`gdalinfo $hgt | grep Size `
gdalwarp -tr 0.001 0.001 -r bilinear $hgt $hgt.ok.tif
rm $hgt
gdal_translate -of GTiff -co COMPRESS=DEFLATE $hgt.ok.tif $hgt >/dev/null 2>/dev/null
rm $hgt.ok.tif

fix_geocoding_frame.sh $frame


