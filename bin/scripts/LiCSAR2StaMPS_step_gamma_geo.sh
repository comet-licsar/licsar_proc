#!/bin/bash
###################################
# step_gamma_geo
# run this script to generate lon.raw and lat.raw
# 
# Dependencies: GMT and GAMMA software
#
# Author: Pablo J. Gonzalez
# Date: 12/Nov/2015
# Last modification: 22/01/2016 Write into geo/ directory (compatibility with gamma2stamps)
###################################

echo " "
PRG=`basename "$0"`
VER="StaMPS"
AUT="Pablo J. Gonzalez (p.j.gonzalez[at]leeds.ac.uk)"
echo "$PRG, $VER, $AUT"
echo " "

if [ -z "$1" ]; then
cat << __EOFHD

  DOCUMENTATION:

          $PRG dem_par lookuptable width length

      The program generates lon.raw and lat.raw files from
      an existing differential interferogram stack made with 
      GAMMA processing software.

            dem_par, dem parameter file of cropped DEM 
        lookuptable, look-up table from geo to radar coordinates
              width, width of the interferograms
             length, length of the interferograms

  EXAMPLE:

      $PRG geo/EQA.dem_par geo/20150919.lt_fine 13458 13579

  DEPENDENCIES:

      awk, GMT [grdmath and grd2xyz], GAMMA [swap_bytes and geocode] http://www.gamma-rs.ch/

__EOFHD
  exit 1
fi

dem_par_file=$1        # e.g., geo/EQA.dem_par
input_lookuptable=$2   # e.g., geo/$masterdate.lt_fine
width=$3               # Width of the interferogram
length=$4              # Length of the interferogram

###########################################################################
# Swap back for Andy
master_date=`ls ../geo/*.hgt | awk '{if(NR==1) print substr($1,5,8)}'`

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
grdmath -R${lon}/${lon1}/${lat1}/${lat} -I${width_dem}+/${length_dem}+ X = geo.grd
grd2xyz geo.grd -ZTLf > geo.raw
swap_bytes geo.raw geolon.raw 4
#rashgt geolon.raw - $width_dem - - - 5 5 0.05 - - - lon_dem.ras # To check results
echo "geocode ${input_lookuptable} geolon.raw ${width_dem} geo/${master_date}.lon ${width} ${length} 2 0 "
geocode ${input_lookuptable} geolon.raw ${width_dem} geo/${master_date}.lon ${width} ${length} 2 0 
swap_bytes geo/${master_date}.lon geo/${master_date}.lon.raw 4
#rashgt tmp.raw - $width - - - 5 5 0.05 - - - lon.ras # To check results

# Latitude file
grdmath -R${lon}/${lon1}/${lat1}/${lat} -I${width_dem}+/${length_dem}+ Y = geo.grd  
grd2xyz geo.grd -ZTLf > geo.raw
swap_bytes geo.raw geolat.raw 4
echo "geocode ${input_lookuptable} geolat.raw ${width_dem} geo/${master_date}.lat ${width} ${length} 2 0 "
geocode ${input_lookuptable} geolat.raw ${width_dem} geo/${master_date}.lat ${width} ${length} 2 0 
swap_bytes geo/${master_date}.lat geo/${master_date}.lat.raw 4
#rashgt tmp.raw - $width - - - 5 5 0.05 - - - lat.ras # To check results

# Cleaning
rm -f geo.raw geolon.raw geolat.raw geo.grd
