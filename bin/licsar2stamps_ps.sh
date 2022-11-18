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
  
__EOFHD
  exit 1
fi
#maxdate=20180501

source $LiCSARpath/lib/LiCSAR_bash_lib.sh

master=`get_master`
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
geocode ${input_lookuptable} geolon.raw ${width_dem} geo/${master_date}.lon ${width} ${length} 2 0 >/dev/null
#swap_bytes geo/${master_date}.lon geo/${master_date}.lon.raw 4 >/dev/null
#rashgt tmp.raw - $width - - - 5 5 0.05 - - - lon.ras # To check results

# Latitude file
gmt grdmath -R${lon}/${lon1}/${lat1}/${lat} -I${width_dem}+/${length_dem}+ Y = geo.grd  
gmt grd2xyz geo.grd -ZTLf > geo.raw
swap_bytes geo.raw geolat.raw 4 >/dev/null 
#rashgt geolat.raw - $width_dem - - - 5 5 0.05 - - - lat_dem.ras # To check results
echo "geocode ${input_lookuptable} geolat.raw ${width_dem} geo/${master_date}.lat ${width} ${length} 2 0 "
geocode ${input_lookuptable} geolat.raw ${width_dem} geo/${master_date}.lat ${width} ${length} 2 0 >/dev/null 
#swap_bytes geo/${master_date}.lat geo/${master_date}.lat.raw 4 >/dev/null 
#rashgt tmp.raw - $width - - - 5 5 0.05 - - - lat.ras # To check results


# stamps gamma expects big endian (?)

mkdir -p INSAR_$master/geo
mkdir INSAR_$master/rslc INSAR_$master/diff0
mv geo/${master_date}.lat INSAR_$master/geo/lat
mv geo/${master_date}.lon INSAR_$master/geo/lon

# make 1-1 hgt
geocode ${input_lookuptable} tostamps.demseg ${width_dem} INSAR_$master/geo/${master}_dem.rdc ${width} ${length} 1 0 >/dev/null 
hgt=INSAR_$master/geo/${master}_dem.rdc




#cp geo/${master}.hgt INSAR_$master/geo/${master}_dem.rdc
#cp geo/${master}.diff_par INSAR_$master/geo/${master}.diff_par

echo "copying rslcs"
for x in `ls RSLC/*/*.rslc`; do cp $x $x.par INSAR_$master/rslc/.; done


echo "generating/copying ifgs and their baselines"
if [ ! -z $maxdate ]; then
for x in `ls RSLC`; do 
 if [ $x != $master ] && [ $x -lt $maxdate ]; then
  echo $x >> slist.txt
 fi
done
else
 ls RSLC >slist.txt
 sed -i '/'$master'/d' slist.txt
fi

for x in `cat slist.txt`; do
 pair=$master'_'$x
 echo 'generating '$pair
 #ifgdir=IFG/$pair
 ifg=INSAR_$master/diff0/$pair.diff
 #create_diff_par RSLC/$master/$master.rslc.par RSLC/$x/$x.rslc.par $ifg'_diffpar'
 create_offset RSLC/$master/$master.rslc.par RSLC/$x/$x.rslc.par $pair.off 1 1 1 0 >/dev/null # for PS - keep no multilooking
 phase_sim_orb RSLC/$master/$master.rslc.par RSLC/$x/$x.rslc.par $pair.off $hgt $pair.simorb >/dev/null
 #SLC_diff_intf RSLC/$master/$master.rslc RSLC/$x/$x.rslc RSLC/$master/$master.rslc.par RSLC/$x/$x.rslc.par $pair.off $pair.simorb $ifg 1 1 0 0 >/dev/null
 SLC_intf2 RSLC/$master/$master.rslc RSLC/$x/$x.rslc RSLC/$master/$master.rslc.par RSLC/$x/$x.rslc.par - - - - $ifg.2 INSAR_$master/diff0/$pair.cc 1 1 - - - - $pair.simorb
 #cc_wave $ifg
 base_init RSLC/$master/$master.rslc.par RSLC/$x/$x.rslc.par $pair.off $ifg INSAR_$master/diff0/$pair.base 0 >/dev/null
 rm $pair.off $pair.simorb
done

# Cleaning
#rm -f geo.raw geolon.raw geolat.raw geo.grd
