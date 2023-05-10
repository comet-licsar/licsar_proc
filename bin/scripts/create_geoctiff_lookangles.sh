#!/bin/bash

procdir=$1
master=$2

if [ ! -d "${procdir}/GEOC" ]; then mkdir ${procdir}/GEOC; fi
if [ ! -d "${procdir}/GEOC/lookangles" ]; then mkdir ${procdir}/GEOC/lookangles; fi

lat=`awk '$1 == "corner_lat:" {print $2}' ${procdir}/geo/EQA.dem_par`
lon=`awk '$1 == "corner_lon:" {print $2}' ${procdir}/geo/EQA.dem_par`
latstep=`awk '$1 == "post_lat:" {print $2}' ${procdir}/geo/EQA.dem_par`
lonstep=`awk '$1 == "post_lon:" {print $2}' ${procdir}/geo/EQA.dem_par`
length_dem=`awk '$1 == "nlines:" {print $2}' ${procdir}/geo/EQA.dem_par`
width_dem=`awk '$1 == "width:" {print $2}' ${procdir}/geo/EQA.dem_par`
reducfac_dem=`echo $width_dem | awk '{if(int($1/2000) > 1) print int($1/2000); else print 1}'`
width=`awk '$1 == "range_samples:" {print $2}' ${procdir}/SLC/$master/$master.slc.mli.par`; #echo $width
lat1=`echo $lat $latstep $length_dem | awk '{printf "%7f", ($1+($2*($3-1)));}'` # Substract one because width starts at zero
lon1=`echo $lon $lonstep $width_dem | awk '{printf "%7f", ($1+($2*($3-1)));}'`
latstep=`awk '$1 == "post_lat:" {if($2<0) printf "%7f", -1*$2; else printf "%7f", $2}' ${procdir}/geo/EQA.dem_par`
lonstep=`awk '$1 == "post_lon:" {if($2<0) printf "%7f", -1*$2; else printf "%7f", $2}' ${procdir}/geo/EQA.dem_par`
# Because no wavelength is reported in master.rmli.par file, we calculated here according to the radar frequency (IN CENTIMETERS)
# Frequency = (C / Wavelength), Where: Frequency: Frequency of the wave in hertz (hz). C: Speed of light (29,979,245,800 cm/sec (3 x 10^10 approx))
lambda=`awk '$1 == "radar_frequency:" {print 29979245800/$2}' ${procdir}/SLC/$master/$master.slc.mli.par`;

echo " Running doGeocoding step "
echo "   check doGeocoding.log if something goes wrong "
logfile=${procdir}/13_doGeocoding.log
rm -f $logfile

echo "   Geocoding results for lookangles." ;
#psi and incidence
if [ ! -f ${procdir}/geo/$master.lt_fine ]; then
 ln -s ${procdir}/geo/$master.lt ${procdir}/geo/$master.lt_fine
fi

if [ -e ${procdir}/geo/psi ]; then
  # Geocode the data
  geocode_back ${procdir}/geo/psi $width_dem ${procdir}/geo/$master.lt_fine ${procdir}/GEOC/lookangles/$master.geo.psi ${width_dem} ${length_dem} 1 0 >> $logfile
  # Convert to geotiff
  #data2geotiff ${procdir}/geo/EQA.dem_par ${procdir}/GEOC/lookangles/$master.geo.psi 2 ${procdir}/GEOC/lookangles/$master.geo.psi.tif 0.0
fi
if [ -e ${procdir}/geo/inc ]; then
  # Geocode the data
  geocode_back ${procdir}/geo/inc $width_dem ${procdir}/geo/$master.lt_fine ${procdir}/GEOC/lookangles/$master.geo.inc ${width_dem} ${length_dem} 1 0 >> $logfile
  # Convert to geotiff
  #data2geotiff ${procdir}/geo/EQA.dem_par ${procdir}/GEOC/lookangles/$master.geo.inc 2 ${procdir}/GEOC/lookangles/$master.geo.inc.tif 0.0
fi
#E-N-U
if [ -e ${procdir}/geo/E ]; then
  # Geocode the data
  geocode_back ${procdir}/geo/E $width ${procdir}/geo/$master.lt_fine ${procdir}/GEOC/lookangles/$master.geo.E ${width_dem} ${length_dem} 1 0 >> $logfile
  # Convert to geotiff
  data2geotiff ${procdir}/geo/EQA.dem_par ${procdir}/GEOC/lookangles/$master.geo.E 2 ${procdir}/GEOC/lookangles/$master.geo.E.tif 0.0
fi
if [ -e ${procdir}/geo/N ]; then
  # Geocode the data
  geocode_back ${procdir}/geo/N $width ${procdir}/geo/$master.lt_fine ${procdir}/GEOC/lookangles/$master.geo.N ${width_dem} ${length_dem} 1 0 >> $logfile
  # Convert to geotiff
  data2geotiff ${procdir}/geo/EQA.dem_par ${procdir}/GEOC/lookangles/$master.geo.N 2 ${procdir}/GEOC/lookangles/$master.geo.N.tif 0.0
fi
if [ -e ${procdir}/geo/U ]; then
  # Geocode the data
  geocode_back ${procdir}/geo/U $width ${procdir}/geo/$master.lt_fine ${procdir}/GEOC/lookangles/$master.geo.U ${width_dem} ${length_dem} 1 0 >> $logfile
  # Convert to geotiff
  data2geotiff ${procdir}/geo/EQA.dem_par ${procdir}/GEOC/lookangles/$master.geo.U 2 ${procdir}/GEOC/lookangles/$master.geo.U.tif 0.0
fi

#hgt - not a 'look angle', but useful for LiCSBAS etc
if [ -e ${procdir}/geo/$master.hgt ]; then
 geocode_back ${procdir}/geo/$master.hgt $width ${procdir}/geo/$master.lt_fine ${procdir}/GEOC/lookangles/$master.geo.hgt ${width_dem} ${length_dem} 1 0 >> $logfile
 data2geotiff ${procdir}/geo/EQA.dem_par ${procdir}/GEOC/lookangles/$master.geo.hgt 2 ${procdir}/GEOC/lookangles/$master.geo.hgt.tif 0.0
fi
