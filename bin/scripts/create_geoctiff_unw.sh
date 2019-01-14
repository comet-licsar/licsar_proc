#!/bin/bash

procdir=$1
master=$2
ifg=$3

echo "Processing dir: $procdir"
echo "Master image: $master"
echo "Interferogram: $ifg"

width=`awk '$1 == "range_samples:" {print $2}' ${procdir}/RSLC/$master/$master.rslc.mli.par`; #echo $width
length=`awk '$1 == "azimuth_lines:" {print $2}' ${procdir}/RSLC/$master/$master.rslc.mli.par`; #echo $length
reducfac=`echo $width | awk '{if(int($1/16000) > 1) print int($1/16000); else print 1}'`

#if [ -d "${procdir}/GEOC" ]; then rm -rf ${procdir}/GEOC; fi
if [ ! -d "${procdir}/GEOC" ]; then mkdir ${procdir}/GEOC; fi

lat=`awk '$1 == "corner_lat:" {print $2}' ${procdir}/geo/EQA.dem_par`
lon=`awk '$1 == "corner_lon:" {print $2}' ${procdir}/geo/EQA.dem_par`
latstep=`awk '$1 == "post_lat:" {print $2}' ${procdir}/geo/EQA.dem_par`
lonstep=`awk '$1 == "post_lon:" {print $2}' ${procdir}/geo/EQA.dem_par`
length_dem=`awk '$1 == "nlines:" {print $2}' ${procdir}/geo/EQA.dem_par`
width_dem=`awk '$1 == "width:" {print $2}' ${procdir}/geo/EQA.dem_par`
reducfac_dem=`echo $width_dem | awk '{if(int($1/2000) > 1) print int($1/2000); else print 1}'`
lat1=`echo $lat $latstep $length_dem | awk '{printf "%7f", ($1+($2*($3-1)));}'` # Substract one because width starts at zero
lon1=`echo $lon $lonstep $width_dem | awk '{printf "%7f", ($1+($2*($3-1)));}'`
latstep=`awk '$1 == "post_lat:" {if($2<0) printf "%7f", -1*$2; else printf "%7f", $2}' ${procdir}/geo/EQA.dem_par`
lonstep=`awk '$1 == "post_lon:" {if($2<0) printf "%7f", -1*$2; else printf "%7f", $2}' ${procdir}/geo/EQA.dem_par`
# Because no wavelength is reported in master.rmli.par file, we calculated here according to the radar frequency (IN CENTIMETERS)
# Frequency = (C / Wavelength), Where: Frequency: Frequency of the wave in hertz (hz). C: Speed of light (29,979,245,800 cm/sec (3 x 10^10 approx))
lambda=`awk '$1 == "radar_frequency:" {print 29979245800/$2}' ${procdir}/RSLC/$master/$master.rslc.mli.par`;

echo " Running doGeocoding step for unwrapped"
echo "   check doGeocoding.log if something goes wrong "
logfile=${procdir}/13_doGeocoding.log
#rm -f $logfile

#mdate=`echo $ifg | awk '{print $2}'`;
#sdate=`echo $ifg | awk '{print $3}'`;
#if [ -d "${procdir}/GEOC/${ifg}" ]; then rm -rf ${procdir}/GEOC/${ifg}; fi
if [ ! -d "${procdir}/GEOC/${ifg}" ]; then mkdir ${procdir}/GEOC/${ifg}; fi
echo "   Geocoding results for inteferogram: ${ifg}" ;
# Unwrapped interferogram
if [ -e ${procdir}/IFG/${ifg}/${ifg}.unw ]; then
  # Replace nan to zeros
  replace_values ${procdir}/IFG/${ifg}/${ifg}.unw nan 0 ${procdir}/IFG/${ifg}/${ifg}.unw0 $width
  mv ${procdir}/IFG/${ifg}/${ifg}.unw0 ${procdir}/IFG/${ifg}/${ifg}.unw
  # Geocode all the data
  geocode_back ${procdir}/IFG/${ifg}/${ifg}.unw $width ${procdir}/geo/$master.lt_fine ${procdir}/GEOC/${ifg}/${ifg}.geo.unw ${width_dem} ${length_dem} 0 0 >> $logfile
  # Convert to geotiff
  data2geotiff ${procdir}/geo/EQA.dem_par ${procdir}/GEOC/${ifg}/${ifg}.geo.unw 2 ${procdir}/GEOC/${ifg}/${ifg}.geo.unw.tif 0.0
  # Create bmps
  #ras_linear ${procdir}/GEOC/${ifg}/${ifg}.geo.unw ${width_dem} - - $reducfac_dem $reducfac_dem - - - ${procdir}/GEOC/${ifg}/${ifg}.geo.unw.bmp >> $logfile
  rasrmg ${procdir}/GEOC/${ifg}/${ifg}.geo.unw ${procdir}/geo/EQA.${master}.slc.mli ${width_dem} - - - $reducfac_dem $reducfac_dem - - - - - ${procdir}/GEOC/${ifg}/${ifg}.geo.unw_blk.bmp >> $logfile
  convert -transparent black -resize 30% ${procdir}/GEOC/${ifg}/${ifg}.geo.unw_blk.bmp ${procdir}/GEOC/${ifg}/${ifg}.geo.unw.bmp

fi


