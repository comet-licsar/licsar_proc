#!/bin/bash

# 2023 

if [ -z $1 ]; then 
 echo "USAGE: provide pair, and keep being in the frame folder"
 echo "e.g. licsar_offset_tracking_pair.sh 20230115_20230127 ... to do px offset tracking between those dates"
 exit
fi

source $LiCSARpath/lib/LiCSAR_bash_lib.sh

m=`echo $1 | cut -d '_' -f1`
s=`echo $1 | cut -d '_' -f2`
pair=$1
master=`get_master`
outdir=IFG/$pair
mkdir -p $outdir
geopairdir=GEOC/$pair
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
create_offset $mpar $spar $outdir/tracking.off 1 1 1 0 >/dev/null

#offset_pwr_tracking $mslc $sslc $mpar $spar tracking.off tracking.offsets tracking.corr 12 4 - 2 -
# avoid oversample... impossible as this is integer, so keep only 2x oversample
#offset_pwr_tracking $mslc $sslc $mpar $spar tracking.off tracking.offsets tracking.corr 32 8 - 1 - >/dev/null
# do 4x oversample
offset_pwr_tracking $mslc $sslc $mpar $spar $outdir/tracking.off $outdir/tracking.offsets $outdir/tracking.corr 64 16 - 2 - >/dev/null
# keeping result in slant range/azimuth - useful for e.g. support in unwrapping:
offset_tracking $outdir/tracking.offsets $outdir/tracking.corr $mpar $outdir/tracking.off $outdir/disp_map $outdir/disp_val 1 - 0 >/dev/null
# now result will be in ground range/azimuth:
#offset_tracking tracking.offsets tracking.corr $mpar tracking.off disp_map disp_val 2 - 0 >/dev/null
widthoff=`grep range_samples $outdir/tracking.off | awk '{print $2}'`
lenoff=`grep azimuth_samples $outdir/tracking.off | awk '{print $2}'`
# extract both range and azi displacements
cpx_to_real $outdir/disp_map $outdir/disp_map.rng $widthoff 0 >/dev/null
cpx_to_real $outdir/disp_map $outdir/disp_map.azi $widthoff 1 >/dev/null
# resample towards orig size
python3 -c "import cv2; import numpy as np; a = np.fromfile('"$outdir/disp_map.rng"', dtype=np.float32).byteswap().reshape(("$lenoff","$widthoff")); cv2.resize(a,dsize=("$mliwid","$mlilen"), interpolation=cv2.INTER_LINEAR).byteswap().tofile('"$outdir/$pair.rng"')" 
python3 -c "import cv2; import numpy as np; a = np.fromfile('"$outdir/disp_map.azi"', dtype=np.float32).byteswap().reshape(("$lenoff","$widthoff")); cv2.resize(a,dsize=("$mliwid","$mlilen"), interpolation=cv2.INTER_LINEAR).byteswap().tofile('"$outdir/$pair.azi"')" 

#now geocode it
mkdir -p $geopairdir
dempar=geo/EQA.dem_par
demwid=`grep width $dempar | gawk {'print $2'}`
geolt=geo/$master.lt_fine
geocode_back $outdir/$pair.rng $mliwid $geolt $geopairdir/$pair.rng.geo $demwid
geocode_back $outdir/$pair.azi $mliwid $geolt $geopairdir/$pair.azi.geo $demwid
#geocode_back disp_map.rng $widthoff $geolt disp_map.rng.geo $demwid
data2geotiff $dempar $geopairdir/$pair.rng.geo 2 $geopairdir/$pair.rng.geo.tif
data2geotiff $dempar $geopairdir/$pair.azi.geo 2 $geopairdir/$pair.azi.geo.tif
#data2geotiff $dempar disp_map.rng.geo 2 $s.rng.geo2.tif
chmod 777 $geopairdir/$pair.rng.geo.tif $geopairdir/$pair.azi.geo.tif

# create previews for the offset geotiffs
