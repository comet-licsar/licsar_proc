#!/bin/bash

# 2023 

if [ -z $1 ]; then 
 echo "USAGE: provide pair, and keep being in the frame folder"
 echo "e.g. licsar_offset_tracking_pair.sh 20230115_20230127 ... to do px offset tracking between those dates"
 exit
fi

unset GOMP_CPU_AFFINITY KMP_AFFINITY   # might help with some parallelisation scaling (?)

source $LiCSARpath/lib/LiCSAR_bash_lib.sh

m=`echo $1 | cut -d '_' -f1`
s=`echo $1 | cut -d '_' -f2`
pair=$1
master=`get_master`
frame=`pwd`
frame=`basename $frame`
outdir=IFG/$pair
mkdir -p $outdir
geopairdir=GEOC/$pair

if [ -f $geopairdir/$pair.geo.rng.tif ]; then
  echo "the range offsets already exist, cancelling"
  exit
fi

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
echo "performing pixel offset tracking"
date
#offset_pwr_tracking $mslc $sslc $mpar $spar $outdir/tracking.off $outdir/tracking.offsets $outdir/tracking.corr 64 16 - 2 - >/dev/null
#time offset_pwr_tracking $mslc $sslc $mpar $spar $outdir/tracking.off $outdir/tracking.offsets $outdir/tracking.corr 128 32 - 2 - #>/dev/null
# deramp first??
if [ ! -f tab/$master'_tab' ]; then
  createSLCtab SLC/$master/$master slc 1 3 > tab/$master'_tab'
fi
for x in $m $s; do
 if [ ! -f tab/$x'R_tab' ]; then
  createSLCtab RSLC/$x/$x rslc 1 3 > tab/$x'R_tab'
 fi
 extd='.deramp'
 if [ ! -f RSLC/$x/$x.rslc.deramp ]; then
  if [ -z `which ScanSAR_deramp_2nd.py 2>/dev/null` ]; then echo "WARNING, old GAMMA - no deramping (but ok if no oversampling)";
   extd=''
  else
   echo "deramping "$x". ETA: 1 minute"
   ScanSAR_deramp_2nd.py tab/$x'R_tab' $x tab/$master'_tab' 20 4 1 >/dev/null
   mv $x.rslc.deramp $x.rslc.deramp.par RSLC/$x/. 
  fi
 fi
done

# only 1 oversample
time offset_pwr_tracking RSLC/$m/$m.rslc$extd RSLC/$s/$s.rslc$extd RSLC/$m/$m.rslc$extd.par RSLC/$s/$s.rslc$extd.par $outdir/tracking.off $outdir/tracking.offsets $outdir/tracking.corr 128 64 - 2 0.1 40 16 - - - - - - 0 1 - - $outdir/tracking.corrstd >/dev/null
# 2^2 oversample
#time offset_pwr_tracking RSLC/$m/$m.rslc.deramp RSLC/$s/$s.rslc.deramp RSLC/$m/$m.rslc.deramp.par RSLC/$s/$s.rslc.deramp.par $outdir/tracking.off $outdir/tracking.offsets $outdir/tracking.corr 64 32 - 2 - >/dev/null
# after Yasser's check: actually gives very very similar result as without deramping! but it is correct to deramp - as only then we can properly oversample, as i tested.
# so will keep deramp, but only 2^1 oversample, since the 2^2 took much more time and no visible improvement, except for higher resolution
#time offset_pwr_tracking $mslc $sslc $mpar $spar $outdir/tracking.off $outdir/tracking.offsets $outdir/tracking.corr 64 32 - 1 - #>/dev/null

echo "done: "
date
echo "extracting the data"
date
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
python3 -c "import cv2; import numpy as np; a = np.fromfile('"$outdir/tracking.corr"', dtype=np.float32).byteswap().reshape(("$lenoff","$widthoff")); cv2.resize(a,dsize=("$mliwid","$mlilen"), interpolation=cv2.INTER_LINEAR).byteswap().tofile('"$outdir/$pair.offsettracking.corr"')" 
python3 -c "import cv2; import numpy as np; a = np.fromfile('"$outdir/tracking.corrstd"', dtype=np.float32).byteswap().reshape(("$lenoff","$widthoff")); cv2.resize(a,dsize=("$mliwid","$mlilen"), interpolation=cv2.INTER_LINEAR).byteswap().tofile('"$outdir/$pair.offsettracking.corrstd"')" 

date
#now geocode it
mkdir -p $geopairdir
dempar=geo/EQA.dem_par
demwid=`grep width $dempar | gawk {'print $2'}`
geolt=geo/$master.lt_fine
echo "geocoding"
geocode_back $outdir/$pair.rng $mliwid $geolt $geopairdir/$pair.rng.geo $demwid >/dev/null
geocode_back $outdir/$pair.azi $mliwid $geolt $geopairdir/$pair.azi.geo $demwid >/dev/null
geocode_back $outdir/$pair.offsettracking.corr $mliwid $geolt $geopairdir/$pair.tracking_corr.geo $demwid >/dev/null
geocode_back $outdir/$pair.offsettracking.corrstd $mliwid $geolt $geopairdir/$pair.tracking_corrstd.geo $demwid >/dev/null

#geocode_back disp_map.rng $widthoff $geolt disp_map.rng.geo $demwid
data2geotiff $dempar $geopairdir/$pair.rng.geo 2 $geopairdir/$pair.geo.rng.tif >/dev/null
data2geotiff $dempar $geopairdir/$pair.azi.geo 2 $geopairdir/$pair.geo.azi.tif >/dev/null
data2geotiff $dempar $geopairdir/$pair.tracking_corr.geo 2 $geopairdir/$pair.geo.tracking_corr.tif >/dev/null
data2geotiff $dempar $geopairdir/$pair.tracking_corrstd.geo 2 $geopairdir/$pair.geo.tracking_corrstd.tif >/dev/null

#data2geotiff $dempar disp_map.rng.geo 2 $s.rng.geo2.tif
chmod 777 $geopairdir/$pair.geo.rng.tif $geopairdir/$pair.geo.azi.tif $geopairdir/$pair.geo.tracking_cor*.tif
rm $geopairdir/$pair.rng.geo $geopairdir/$pair.azi.geo $geopairdir/$pair.tracking_corr.geo $geopairdir/$pair.tracking_corrstd.geo

# create previews for the offset geotiffs
create_preview_offsets $geopairdir/$pair.geo.rng.tif $frame 10
create_preview_offsets $geopairdir/$pair.geo.azi.tif $frame 10
