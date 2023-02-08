#!/bin/bash

# 2023 

if [ -z $1 ]; then 
 echo "USAGE: provide pair, and keep being in the frame folder"
 echo "e.g. licsar_offset_tracking_pair.sh 20230115_20230127 ... to do px offset tracking between those dates"
 exit
fi

source $LiCSARpath/lib/LiCSAR_bash_lib.sh

m=`echo $1 | cut '_' -f1`
s=`echo $1 | cut '_' -f2`
master=`get_master`

# idea is - perform the offset tracking without oversampling, with rg=32, az=8
# then resample to the MLI dimensions (20/4 multilooking)
# and geocode to geotiff
rslcdir=`pwd`/RSLC
mmli=$rslcdir/$m/$m.rslc.mli
mpar=$rslcdir/$m/$m.rslc.par
mslc=$rslcdir/$m/$m.rslc
mliwid=`grep range_samples $mmli.par | gawk {'print $2'}`
mlilen=`grep azimuth_lines $mmli.par | gawk {'print $2'}`
spar=$rslcdir/$s/$s.rslc.par
sslc=$rslcdir/$s/$s.rslc
create_offset $mpar $spar tracking.off 1 1 1 0 >/dev/null

#offset_pwr_tracking $mslc $sslc $mpar $spar tracking.off tracking.offsets tracking.corr 12 4 - 2 -
# avoid oversample... impossible as this is integer, so keep only 2x oversample
offset_pwr_tracking $mslc $sslc $mpar $spar tracking.off tracking.offsets tracking.corr 32 8 - 1 - >/dev/null
offset_tracking tracking.offsets tracking.corr $mpar tracking.off disp_map disp_val 1 - 0 >/dev/null
widthoff=`grep range_samples tracking.off | awk '{print $2}'`
lenoff=`grep azimuth_samples tracking.off | awk '{print $2}'`
cpx_to_real disp_map disp_map.rng $widthoff 0 >/dev/null
# resample towards orig size
python3 -c "import cv2; import numpy as np; a = np.fromfile('disp_map.rng', dtype=np.float32).byteswap().reshape(("$lenoff","$widthoff")); cv2.resize(a,dsize=("$mliwid","$mlilen"), interpolation=cv2.INTER_LINEAR).byteswap().tofile('"$s.rng"')" 

#now geocode it
dempar=geo/EQA.dem_par
demwid=`grep width $dempar | gawk {'print $2'}`
geolt=geo/$master.lt_fine
geocode_back $s.rng $mliwid $geolt $s.rng.geo $demwid
geocode_back disp_map.rng $widthoff $geolt disp_map.rng.geo $demwid
data2geotiff $dempar $s.rng.geo 2 $s.rng.geo.tif
data2geotiff $dempar disp_map.rng.geo 2 $s.rng.geo2.tif

# right, so.. just store it to the GEOC folder..
mkdir -p GEOC/$1
mv $s.rng.geo.tif GEOC/$1/$1.rng.geo.tif
mv $s.rng2.geo.tif GEOC/$1/$1.rng2.geo.tif
