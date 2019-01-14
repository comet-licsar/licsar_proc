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

echo " Running doGeocoding step "
echo "   check doGeocoding.log if something goes wrong "
logfile=${procdir}/13_doGeocoding.log
rm -f $logfile

#mdate=`echo $ifg | awk '{print $2}'`;
#sdate=`echo $ifg | awk '{print $3}'`;
if [ -d "${procdir}/GEOC/${ifg}" ]; then rm -rf ${procdir}/GEOC/${ifg}; fi
if [ ! -d "${procdir}/GEOC/${ifg}" ]; then mkdir ${procdir}/GEOC/${ifg}; fi
echo "   Geocoding results for inteferogram: ${ifg}" ;
# Unfiltered interferogram
if [ -e ${procdir}/IFG/${ifg}/${ifg}.diff ]; then
  # Exctract the mag and phase
  #/home/users/ehatton/insartools/include/cpx2real.py ${procdir}/IFG/${ifg}/${ifg}.diff ${procdir}/IFG/${ifg}/${ifg}.diff_mag $width 3
  #/home/users/ehatton/insartools/include/cpx2real.py ${procdir}/IFG/${ifg}/${ifg}.diff ${procdir}/IFG/${ifg}/${ifg}.diff_pha $width 4
  cpx_to_real ${procdir}/IFG/${ifg}/${ifg}.diff ${procdir}/IFG/${ifg}/${ifg}.diff_mag $width 3
  cpx_to_real ${procdir}/IFG/${ifg}/${ifg}.diff ${procdir}/IFG/${ifg}/${ifg}.diff_pha $width 4
  # Geocode all the data
  geocode_back ${procdir}/IFG/${ifg}/${ifg}.diff $width ${procdir}/geo/$master.lt_fine ${procdir}/GEOC/${ifg}/${ifg}.geo.diff ${width_dem} ${length_dem} 1 1 >> $logfile
  geocode_back ${procdir}/IFG/${ifg}/${ifg}.diff_mag $width ${procdir}/geo/$master.lt_fine ${procdir}/GEOC/${ifg}/${ifg}.geo.diff_mag ${width_dem} ${length_dem} 1 0 >> $logfile
  geocode_back ${procdir}/IFG/${ifg}/${ifg}.diff_pha $width ${procdir}/geo/$master.lt_fine ${procdir}/GEOC/${ifg}/${ifg}.geo.diff_pha ${width_dem} ${length_dem} 0 0 >> $logfile
  # Convert to geotiff
  #data2geotiff ${procdir}/geo/EQA.dem_par ${procdir}/GEOC/${ifg}/${ifg}.geo.diff 4 ${procdir}/GEOC/${ifg}/${ifg}.geo.diff.tif 0.0
  data2geotiff ${procdir}/geo/EQA.dem_par ${procdir}/GEOC/${ifg}/${ifg}.geo.diff_mag 2 ${procdir}/GEOC/${ifg}/${ifg}.geo.diff_mag.tif 0.0
  data2geotiff ${procdir}/geo/EQA.dem_par ${procdir}/GEOC/${ifg}/${ifg}.geo.diff_pha 2 ${procdir}/GEOC/${ifg}/${ifg}.geo.diff_pha.tif 0.0
  # Create bmps
  rasmph_pwr ${procdir}/GEOC/${ifg}/${ifg}.geo.diff ${procdir}/geo/EQA.${master}.slc.mli ${width_dem} - - - $reducfac_dem $reducfac_dem - - - ${procdir}/GEOC/${ifg}/${ifg}.geo.diff_blk.bmp >> $logfile
  convert -transparent black -resize 30% ${procdir}/GEOC/${ifg}/${ifg}.geo.diff_blk.bmp ${procdir}/GEOC/${ifg}/${ifg}.geo.diff.bmp
  #rasmph_pwr ${procdir}/GEOC/${ifg}/${ifg}.geo.diff_mag ${procdir}/geo/EQA.${master}.slc.mli ${width_dem} - - - $reducfac_dem $reducfac_dem - - - ${procdir}/GEOC/${ifg}/${ifg}.geo.diff_mag.bmp >> $logfile
  #rasmph_pwr ${procdir}/GEOC/${ifg}/${ifg}.geo.diff_pha ${procdir}/geo/EQA.${master}.slc.mli ${width_dem} - - - $reducfac_dem $reducfac_dem - - - ${procdir}/GEOC/${ifg}/${ifg}.geo.diff_pha.bmp >> $logfile
  raspwr ${procdir}/GEOC/${ifg}/${ifg}.geo.diff_mag ${width_dem} - - $reducfac_dem $reducfac_dem - - - ${procdir}/GEOC/${ifg}/${ifg}.geo.diff_mag.bmp >> $logfile
  ras_linear ${procdir}/GEOC/${ifg}/${ifg}.geo.diff_pha ${width_dem} - - $reducfac_dem $reducfac_dem - - - ${procdir}/GEOC/${ifg}/${ifg}.geo.diff_pha.bmp >> $logfile
fi

# Unfiltered coherence
if [ -e ${procdir}/IFG/${ifg}/${ifg}.cc ]; then
  # Geocode
  geocode_back ${procdir}/IFG/${ifg}/${ifg}.cc $width ${procdir}/geo/$master.lt_fine ${procdir}/GEOC/${ifg}/${ifg}.geo.cc ${width_dem} ${length_dem} 1 0 >> $logfile
  # Convert to geotiff
  data2geotiff ${procdir}/geo/EQA.dem_par ${procdir}/GEOC/${ifg}/${ifg}.geo.cc 2 ${procdir}/GEOC/${ifg}/${ifg}.geo.cc.tif 0.0
  # create bmps
  rascc ${procdir}/GEOC/${ifg}/${ifg}.geo.cc - ${width_dem} - - - $reducfac_dem $reducfac_dem 0 1 - - - ${procdir}/GEOC/${ifg}/${ifg}.geo.cc_blk.bmp >> $logfile
  # Need to remove the black border, but the command below is no good as it removes the black parts of the coherence!
  #convert -transparent black -resize 30% ${procdir}/GEOC/${ifg}/${ifg}.geo.cc_blk.bmp ${procdir}/GEOC/${ifg}/${ifg}.geo.cc.bmp
  convert -resize 30% ${procdir}/GEOC/${ifg}/${ifg}.geo.cc_blk.bmp ${procdir}/GEOC/${ifg}/${ifg}.geo.cc.bmp
fi

