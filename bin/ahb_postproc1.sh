#!/bin/bash

if [ -z $1 ]; then
 echo "This script will perform extra processing after licsar2licsbas.sh"
 echo "PLEASE run it inside your frame directory (where you have GEOCmlX and TS_GEOCmlX directories)"
 echo "parameters: GEOCDIR frameID"
 echo "e.g. GEOCml10GACOSmask 155D_02611_050400"
 exit
fi

geocd=$1
frame=$2

if [ ! -d TS_$geocd ]; then
  echo "Seems results do not exist here, cancelling"
  exit
fi

if [ ! -f $frame.vstd_scaled.geo.tif ]; then
 LiCSBAS_vel_plate_motion.py -t TS_$geocd -f $frame -o $frame.vel_filt.mskd.eurasia.geo.tif --vstd_fix
 cp TS_$geocd/results/vstd_scaled.tif $frame.vstd_scaled.geo.tif
fi

for x in U E N hgt; do
  if [ ! -f $frame.$x.geo.tif ]; then
    LiCSBAS_flt2geotiff.py -i $geocd/$x -p $geocd/EQA.dem_par -o $frame.$x.geo.tif
  fi
done
