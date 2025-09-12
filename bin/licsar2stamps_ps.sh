#!/bin/bash
###################################
# run this script to prepare dataset for STAMPS PS processing
# 
# Dependencies: GMT and GAMMA software
#
# Author: Milan Lazecky
###################################

if [ -z "$1" ]; then
cat << __EOFHD

  based on Pablo Gonzalez code from 2015-11-12, but after major modification - this script allows exporting frame data for PS processing in stamps
    just run with any one parameter.
  
  testing version. we assume you are in standard LiCSAR processing directory, i.e. with geo, DEM, RSLC folders.

  You can set custom reference epoch, second parameter..

__EOFHD
  exit 1
fi
#maxdate=20180501

source $LiCSARpath/lib/LiCSAR_bash_lib.sh

if [ ! -z $2 ]; then
  master=$2
  if [ ! -d RSLC/$master ]; then echo "no such RSLC: "$master; exit; fi
  echo "setting reference epoch to "$master
else
  master=`get_master`
fi
master_date=$master

# start by making new lt - and skip the 'fine' now

#dem_par_file=geo/EQA.dem_par        # e.g., geo/EQA.dem_par
#input_lookuptable=geo/$master.lt_fine   # e.g., geo/$masterdate.lt_fine
if [ ! -f DEM/dem_crop.dem ]; then
 echo "please include DEM in this folder"
 exit
fi

echo 'preparing data for geocoding'
gc_map RSLC/$master/$master.rslc.par - DEM/dem_crop.dem_par DEM/dem_crop.dem tostamps.demseg.par tostamps.demseg tostamps.lt - - - - - tostamps.inc - - - 0 1 >/dev/null
dem_par_file=tostamps.demseg.par
input_lookuptable=tostamps.lt

width=`grep range_samples RSLC/$master/$master.rslc.par | gawk {'print $2'}`
length=`grep azimuth_lines RSLC/$master/$master.rslc.par | gawk {'print $2'}`
rangeres=`grep range_pixel_spacing RSLC/$master/$master.rslc.par | gawk {'print $2'}`


# Extract geocoding information to generate the lon and lat matrices
lat=`awk '$1 == "corner_lat:" {print $2}' ${dem_par_file}`
lon=`awk '$1 == "corner_lon:" {print $2}' ${dem_par_file}`
latstep=`awk '$1 == "post_lat:" {print $2}' ${dem_par_file}`
lonstep=`awk '$1 == "post_lon:" {print $2}' ${dem_par_file}`
length_dem=`awk '$1 == "nlines:" {print $2}' ${dem_par_file}`
width_dem=`awk '$1 == "width:" {print $2}' ${dem_par_file}`
lat1=`echo $lat $latstep $length_dem | awk '{printf "%7f", ($1+($2*($3-1)));}'` # Substract one because width starts at zero
lon1=`echo $lon $lonstep $width_dem | awk '{printf "%7f", ($1+($2*($3-1)));}'`
latstep=`awk '$1 == "post_lat:" {if($2<0) printf "%7f", -1*$2; else printf "%7f", $2}' ${dem_par_file}`
lonstep=`awk '$1 == "post_lon:" {if($2<0) printf "%7f", -1*$2; else printf "%7f", $2}' ${dem_par_file}`

# From geocoded coordinates to radar coordinates
# Longitude file
gmt grdmath -R${lon}/${lon1}/${lat1}/${lat} -I${width_dem}+/${length_dem}+ X = geo.grd	#make the lon grid
#gmt grdmath -R${lon}/${lon1}/${lat1}/${lat} -I${width}+/${length}+ X = geo.grd	#make the lon grid
gmt grd2xyz geo.grd -ZTLf > geo.raw # take lons
swap_bytes geo.raw geolon.raw 4 >/dev/null # set lons to 4-byte floats
#rashgt geolon.raw - $width_dem - - - 5 5 0.05 - - - lon_dem.ras # To check results
echo "geocode ${input_lookuptable} geolon.raw ${width_dem} geo/${master_date}.lon ${width} ${length} 2 0 "
geocode ${input_lookuptable} geolon.raw ${width_dem} geo/${master_date}.lon ${width} ${length} 0 0 >/dev/null



#swap_bytes geo/${master_date}.lon geo/${master_date}.lon.raw 4 >/dev/null
#rashgt tmp.raw - $width - - - 5 5 0.05 - - - lon.ras # To check results

# Latitude file
gmt grdmath -R${lon}/${lon1}/${lat1}/${lat} -I${width_dem}+/${length_dem}+ Y = geo.grd  
gmt grd2xyz geo.grd -ZTLf > geo.raw
swap_bytes geo.raw geolat.raw 4 >/dev/null 
#rashgt geolat.raw - $width_dem - - - 5 5 0.05 - - - lat_dem.ras # To check results
echo "geocode ${input_lookuptable} geolat.raw ${width_dem} geo/${master_date}.lat ${width} ${length} 2 0 "
geocode ${input_lookuptable} geolat.raw ${width_dem} geo/${master_date}.lat ${width} ${length} 0 0 >/dev/null 
#swap_bytes geo/${master_date}.lat geo/${master_date}.lat.raw 4 >/dev/null 
#rashgt tmp.raw - $width - - - 5 5 0.05 - - - lat.ras # To check results



# let us prepare everything in little endian (also, use adapted mt_prep_gamma - see /nfs/a1/software/StaMPS_v4.1b_ML/bin/mt_prep_gamma in hal)

mkdir -p INSAR_$master/geo
mkdir INSAR_$master/rslc INSAR_$master/diff0
swap_bytes geo/${master_date}.lat INSAR_$master/geo/$master.lat 4 >/dev/null
swap_bytes geo/${master_date}.lon INSAR_$master/geo/$master.lon 4 >/dev/null



cat << EOF >checkthis_correct_lonlat_as_it_might_be_needed.py
# print('If you run this from INSAR* folder, please avoid byteswapping below!')
width=$width
length=$length
lonfile='INSAR_$master/geo/$master.lon'
latfile='INSAR_$master/geo/$master.lat'

import numpy as np
import matplotlib.pyplot as plt
import xarray as xr

lon=np.fromfile(lonfile, dtype=np.float32) #.byteswap()
lat=np.fromfile(latfile, dtype=np.float32) #.byteswap()

lon=lon.reshape((length,width))
lat=lat.reshape((length,width))
lon[lon==0] = np.nan
lat[lat==0] = np.nan

lonxr=xr.DataArray(lon)
#lonxr.plot(); plt.show()
lonxr=lonxr.interpolate_na('dim_0','linear')
lonxr=lonxr.interpolate_na('dim_1','linear')

latxr=xr.DataArray(lat)
latxr=latxr.interpolate_na('dim_0','linear')
latxr=latxr.interpolate_na('dim_1','linear')

lat=latxr.fillna(0).values
lon=lonxr.fillna(0).values

lat.tofile(latfile)
lon.tofile(lonfile)

exit()
EOF

# ok, let's just run this..
echo "updating lon/lat files"
python3 checkthis_correct_lonlat_as_it_might_be_needed.py

# make 1-1 hgt
hgt=geo/${master}_dem.rdc
geocode ${input_lookuptable} tostamps.demseg ${width_dem} $hgt ${width} ${length} 1 0 >/dev/null 
swap_bytes $hgt INSAR_$master/geo/${master}_dem.rdc 4 >/dev/null


#cp geo/${master}.hgt INSAR_$master/geo/${master}_dem.rdc
#cp geo/${master}.diff_par INSAR_$master/geo/${master}.diff_par


echo "selecting dataset"
msize=`ls -al RSLC/$master/$master.rslc | gawk {'print $5'}`
if [ -z $maxdate ]; then maxdate=999999999; fi
for x in `ls RSLC`; do 
 if [ $x != $master ] && [ $x -lt $maxdate ] && [ -f RSLC/$x/$x.rslc ]; then
  ssize=`ls -al RSLC/$x/$x.rslc | gawk {'print $5'}`
  if [ $msize -eq $ssize ]; then
   echo $x >> slist.txt
  fi
 fi
done

###

#master=`ls geo/*rdc | cut -d '/' -f2 | cut -d '_' -f1`
#msize=`ls -al rslc/$master.rslc | gawk {'print $5'}` 
#for xx in `ls rslc/*.rslc`; do
# x=`basename $xx | cut -d '.' -f1`
# ssize=`ls -al rslc/$x.rslc | gawk {'print $5'}` 
# if [ ! $msize -eq $ssize ]; then
#   echo $x
#   rm rslc/$x.rsl* diff0/*$x*
# fi
#done
###

echo "generating/copying ifgs and their baselines, byteswapping rslcs"
mkdir -p stampsifgtemp
for x in `cat slist.txt`; do
 pair=$master'_'$x
 echo 'generating '$pair
 ifg=stampsifgtemp/$pair.diff
 #create_diff_par RSLC/$master/$master.rslc.par RSLC/$x/$x.rslc.par $ifg'_diffpar'
 create_offset RSLC/$master/$master.rslc.par RSLC/$x/$x.rslc.par $pair.off 1 1 1 0 >/dev/null # for PS - keep no multilooking
 phase_sim_orb RSLC/$master/$master.rslc.par RSLC/$x/$x.rslc.par $pair.off $hgt $pair.simorb >/dev/null
 #SLC_diff_intf RSLC/$master/$master.rslc RSLC/$x/$x.rslc RSLC/$master/$master.rslc.par RSLC/$x/$x.rslc.par $pair.off $pair.simorb $ifg 1 1 0 0 >/dev/null
 SLC_intf2 RSLC/$master/$master.rslc RSLC/$x/$x.rslc RSLC/$master/$master.rslc.par RSLC/$x/$x.rslc.par - - - - $ifg stampsifgtemp/$pair.cc 1 1 - - - - $pair.simorb
 #cc_wave $ifg
 base_init RSLC/$master/$master.rslc.par RSLC/$x/$x.rslc.par $pair.off $ifg INSAR_$master/diff0/$pair.base 0 >/dev/null
 if [ `ls -al INSAR_$master/diff0/$pair.base | gawk {'print $5'}` -eq 0 ]; then
  echo "some error with epoch "$x" - skipping"
 else
  # swapping to little endian: ifgs are FCOMPLEX, rslcs are SCOMPLEX, ccs are FLOAT
  swap_bytes $ifg INSAR_$master/diff0/$pair.diff 4 >/dev/null
  swap_bytes stampsifgtemp/$pair.cc INSAR_$master/diff0/$pair.cc 4 >/dev/null
  cp RSLC/$x/$x.rslc.par INSAR_$master/rslc/.
  swap_bytes RSLC/$x/$x.rslc INSAR_$master/rslc/$x.rslc 2 >/dev/null
 fi
 rm $pair.off $pair.simorb
done
# adding also the reference epoch:
cp RSLC/$master/$master.rslc.par INSAR_$master/rslc/.
swap_bytes RSLC/$master/$master.rslc INSAR_$master/rslc/$master.rslc 2 >/dev/null

# cleaning
rm -rf stampsifgtemp



# adding more info
echo 1 > INSAR_$master/slc_osfactor.1.in
echo 0.05546576 > INSAR_$master/lambda.1.in
#echo 2.3295 > rangepixelsize.1.in
echo $rangeres > INSAR_$master/rangepixelsize.1.in

# Cleaning
#rm -f geo.raw geolon.raw geolat.raw geo.grd

pxperpatchsq=500
#width=`cat width.txt`
#master=`ls geo/*lat | cut -d '/' -f2 | cut -d '.' -f1`

let RP=$width/$pxperpatchsq
let AP=$length/$pxperpatchsq
RP=5
AP=1
cat << EOF > INSAR_$master/stamps_proc.sh
source /nfs/a1/software/StaMPS_v4.1b_ML/StaMPS_CONFIG.bash
module load matlab

# mt_prep_gamma $master \`pwd\` 0.4 $RP $AP 100 50 >/dev/null

ls rslc/*.rslc > calamp.in
calamp calamp.in $width calamp.out s 0

# extra check (it helps)
for epoch in \`grep ".rslc 0" calamp.out | gawk {'print \$1'} | rev | cut -d '.' -f2 | cut -d '/' -f1 | rev\`; do
 echo "removing empty "\$epoch
 rm diff0/*$epoch.* rslc/$epoch.*
done

mt_prep_gamma $master \`pwd\` 0.4 $RP $AP 100 50

# matlab -nodesktop -nosplash -r "addpath('/nfs/a1/software/StaMPS_v4.1b_ML/matlab'); addpath('/nfs/a1/software/StaMPS_bjmarfito/matlab'); \
# setparm('ref_radius',300);setparm('ref_centre_lonlat',[18.55752423, 49.84894293]); \

matlab -nodesktop -nosplash -r "addpath('/nfs/a1/software/StaMPS_v4.1b_ML/matlab'); \
getparm; setparm('max_topo_err',50); setparm('gamma_change_convergence',0.01); setparm('gamma_max_iterations',3); \
setparm('weed_time_win',365); setparm('weed_standard_dev',1.3); \
setparm('clap_win',32); setparm('small_baseline_flag','n'); \
setparm('merge_resample_size',0); setparm('small_baseline_flag','n');  \
setparm('unwrap_time',365); setparm('unwrap_gold_alpha',0.8); setparm('unwrap_gold_n_win',16); setparm('unwrap_grid',200); \
stamps(1,6); setparm('unwrap_spatial_cost_func_flag','n'); setparm('subtr_tropo','n'); stamps(7,7); \
ps_calc_ifg_std;a=load('ifgstd2.mat');setparm('scla_drop_i',find(a.ifg_std>50)'); \
stamps(6,6); setparm('unwrap_spatial_cost_func_flag','y'); stamps(7,7); \
setparm('unwrap_spatial_cost_func_flag','n'); stamps(6,7); \
a=load('ifgstd2.mat');setparm('scla_drop_i',find(a.ifg_std>55)'); setparm('unwrap_hold_good_values','y'); \
for u=1:2, u, stamps(6,7); end; \
stamps(6,6); ps_plot('V-do',-1); \
it4s1_stamps2csv; exit"
it4s1_convert2okcsv.sh exported.csv 

EOF
chmod 777 INSAR_$master/stamps_proc.sh
echo "now copy the whole INSAR_ folder to Leeds Uni server and run:"
ls INSAR_$master/stamps_proc.sh

