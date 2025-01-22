#!/bin/bash

# adding flags to control deramping and step sizes of range and azi - Muhammet Nergizci, 2025
# Milan Lazecky, 2023 

if [ -z "$1" ]; then 
 echo "USAGE: provide pair in format YYYYMMDD_YYYYMMDD and run in the frame folder"
 echo "e.g. licsar_offset_tracking_pair.sh 20230115_20230127 (--noderamp) [--rwin <value>] [--awin <value>] [--rstep <value>] [--astep <value>] ... to do px offset tracking between those dates"
 exit
fi

unset GOMP_CPU_AFFINITY KMP_AFFINITY   # might help with some parallelisation scaling (?)

source $LiCSARpath/lib/LiCSAR_bash_lib.sh

# Extract and validate input pair before processing additional arguments
pair=$1
# Extract master and slave dates
m=$(echo "$pair" | cut -d '_' -f1)
s=$(echo "$pair" | cut -d '_' -f2)


# Default values for rwin, awin, rstep, and astep
rwin=128
awin=128
rstep=40
astep=40
novr=2

# Check for the -noderamp flag and other optional arguments
noderamp=false
while [[ $# -gt 0 ]]; do
  case "$1" in
    --noderamp)
      noderamp=true
      echo "Deramping is disabled with --noderamp flag."
      shift
      ;;
    --rwin)
      rwin="$2"
      shift 2
      ;;
    --awin)
      awin="$2"
      shift 2
      ;;
    --rstep)
      rstep="$2"
      shift 2
      ;;
    --astep)
      astep="$2"
      shift 2
      ;;
    --novr)
      novr="$2"
      shift 2
      ;;
    *)
      shift
      ;;
  esac
done

# Get master frame
master=`get_master`
frame=`pwd`
frame=`basename $frame`
if [ `echo $frame | wc -m` != 18 ]; then echo "not in frame folder"; exit; fi

outdir="IFG/$pair"
mkdir -p "$outdir"
geopairdir="GEOC/$pair"

if [ -f "$geopairdir/$pair.geo.rng.tif" ]; then
  echo "The range offsets already exist, cancelling"
  exit
fi


# idea is - perform the offset tracking without oversampling, with rg=32, az=8
# then resample to the MLI dimensions (20/4 multilooking)
# and geocode to geotiff
rslcdir=`pwd`/RSLC
mmli=$rslcdir/$master/$master.rslc.mli
if [ ! -f $mmli ]; then
 echo "additional multilooking..."
 multi_look $rslcdir/$master/$master.rslc $rslcdir/$master/$master.rslc.par $mmli $mmli.par 20 4 >/dev/null 2>/dev/null
fi
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
echo "WARNING, TABS ARE RECREATED BUT MIGHT BE WRONG IF NOT ON 3 SWATHS"
date
#offset_pwr_tracking $mslc $sslc $mpar $spar $outdir/tracking.off $outdir/tracking.offsets $outdir/tracking.corr 64 16 - 2 - >/dev/null
#time offset_pwr_tracking $mslc $sslc $mpar $spar $outdir/tracking.off $outdir/tracking.offsets $outdir/tracking.corr 128 32 - 2 - #>/dev/null
# deramp first??
if [ ! -d tab ]; then mkdir tab; fi
if [ ! -f tab/$master'_tab' ]; then
  createSLCtab_frame SLC/$master/$master slc $frame > tab/$master'_tab'
fi

for x in $m $s; do
  echo "Processing $x..."
  
  # Check and create tab file if it doesn't exist
  if [ ! -f tab/$x'R_tab' ]; then
    createSLCtab_frame RSLC/$x/$x rslc $frame > tab/$x'R_tab'
  fi

  # Default extension
  extd='.deramp'

  # Handle noderamp flag
  if [ "$noderamp" = true ]; then
    extd=''
    echo "Skipping deramping for $x."
  else
    # Perform deramping if not already done
    if [ ! -f RSLC/$x/$x.rslc.deramp ]; then
      if [ -z "$(which ScanSAR_deramp_2nd.py 2>/dev/null)" ]; then 
        echo "WARNING: ScanSAR_deramp_2nd.py not found. Skipping deramping for $x."
        extd=''
      else
        echo "Deramping $x. ETA: 1 minute"
        ScanSAR_deramp_2nd.py tab/$x'R_tab' $x tab/$master'_tab' 20 4 1 >/dev/null
        mv $x.rslc.deramp $x.rslc.deramp.par RSLC/$x/.
        extd='.deramp'
      fi
    fi
  fi
done

echo $extd $rwin $awin $rstep $astep $novr

#offset_pwr_tracking <SLC1> <SLC2> <SLC1_par> <SLC2_par> <OFF_par> <offs> <ccp> ##[rwin] [azwin] [offsets] [n_ovr] [thres] [rstep] [azstep] [rstart] [rstop] [azstart] [azstop] [lanczos] [bw_frac] [deramp] [int_filt] [pflag] [pltflg] [ccs]
time offset_pwr_tracking RSLC/$m/$m.rslc$extd RSLC/$s/$s.rslc$extd RSLC/$m/$m.rslc$extd.par RSLC/$s/$s.rslc$extd.par $outdir/tracking.off $outdir/tracking.offsets $outdir/tracking.corr $rwin $awin - $novr 0.1 $rstep $astep - - - - - - 0 1 - - $outdir/tracking.corrstd >/dev/null
#time offset_pwr_tracking RSLC/$m/$m.rslc$extd RSLC/$s/$s.rslc$extd RSLC/$m/$m.rslc$extd.par RSLC/$s/$s.rslc$extd.par $outdir/tracking.off $outdir/tracking.offsets $outdir/tracking.corr 128 128 - 2 0.1 40 40 - - - - - - 0 1 - - $outdir/tracking.corrstd >/dev/null
#time offset_pwr_tracking RSLC/$m/$m.rslc$extd RSLC/$s/$s.rslc$extd RSLC/$m/$m.rslc$extd.par RSLC/$s/$s.rslc$extd.par $outdir/tracking.off $outdir/tracking.offsets $outdir/tracking.corr 64 32 - 1 0.1 20 10 - - - - - - 0 1 - - $outdir/tracking.corrstd >/dev/null

##Milan's notes:

# ok ok, so now i do only the deramped...
#e.g., for the test_EQ_tur:
#cd /gws/nopw/j04/nceo_geohazards_vol1/projects/LiCS/proc/current/subsets/test_tur_rs/021D
#time offset_pwr_tracking RSLC/$m/$m.rslc RSLC/$s/$s.rslc RSLC/$m/$m.rslc.par RSLC/$s/$s.rslc.par $outdir/tracking.off $outdir/tracking.offsets $outdir/tracking.corr 128 64 - 2 0.2 20 4 #- - - - - - 0 1 - - $outdir/tracking.corrstd

# only 1 oversample
#time offset_pwr_tracking RSLC/$m/$m.rslc$extd RSLC/$s/$s.rslc$extd RSLC/$m/$m.rslc$extd.par RSLC/$s/$s.rslc$extd.par $outdir/tracking.off $outdir/tracking.offsets $outdir/tracking.corr 128 128 - 2 0.1 40 40 - - - - - - 0 1 - - $outdir/tracking.corrstd >/dev/null
#time offset_pwr_tracking RSLC/$m/$m.rslc$extd RSLC/$s/$s.rslc$extd RSLC/$m/$m.rslc$extd.par RSLC/$s/$s.rslc$extd.par $outdir/tracking.off $outdir/tracking.offsets $outdir/tracking.corr 128 64 - 2 0.1 40 16 - - - - - - 0 1 - - $outdir/tracking.corrstd >/dev/null

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
