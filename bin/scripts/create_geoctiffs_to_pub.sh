#!/bin/bash
if [ -z $2 ]; then 
 echo "inputs are: procdir ifg, e.g. \`pwd\` 20160101_20160202";
 echo "optional parameters:"
 echo "-u .... geocode also unfiltered wrapped interferogram"
 echo "-F .... do full resolution previews (needed for KML)"
 echo "-L .... do low resolution geotiffs (500x500 m)"
 echo "-H .... do full resolution geotiffs (50x50 m)"
 echo "some more parameters:"
 echo "-U .... geocode only unwrapped interferogram"
 echo "-I .... geocode only wrapped interferogram"
 echo "-C .... geocode only coherence"
 echo "-M .... geocode only MLI"
 exit;
fi
#what is needed here (can be changed probably) is *.rslc.mli.par of master image !!!

#module load doris
#module load LiCSAR/dev
RESIZE=30
#to have same resolution in PNG files, use:
#RESIZE=100
UNWDO=1
IFGDO=1
COHDO=1
MLIDO=1
HIRES=0
LORES=0
UNFILT=0

while getopts ":uUCFHIML" option; do
 case "${option}" in
  u) UNFILT=1; echo "unfiltered ifg will also be geocoded";
     ;;
  U) IFGDO=0; COHDO=0; MLIDO=0; echo "you do ONLY unw geo now..";
     ;;
  C) IFGDO=0; UNWDO=0; MLIDO=0; echo "you do ONLY coh geo now..";
     ;;
  I) COHDO=0; UNWDO=0; MLIDO=0; echo "you do ONLY ifg geo now..";
     ;;
  M) COHDO=0; UNWDO=0; IFGDO=0; echo "you do ONLY MLIs geo now..";
     ;;
  F) RESIZE=100; echo "previews will be generated in full resolution";
     ;;
  H) HIRES=1; echo "geotiffs will be generated in full resolution (50x50 m)";
     ;;
  L) LORES=1; echo "geotiffs will be generated in low resolution (500x500 m)";
     ;;
 esac
done
shift $((OPTIND -1))


procdir=$1
#master=$2
ifg=$2
geodir=geo

GEOCDIR=GEOC

if [ $HIRES == 1 ]; then
 if [ ! -d $procdir/geo_50m ]; then
  echo "generating high resolution geo files"
  thisdir=`pwd`
  cd $procdir
  submit_geo_hires.py
  geodir=geo_50m
 fi
 if [ -f $procdir/geo_50m/locked ]; then
  echo "seems like the hires geocoding is locked by another process. please wait for it to finish first"
  exit
 fi

 GEOCDIR=GEOC_50m
 mkdir -p $GEOCDIR
 cd $thisdir
fi

if [ $LORES == 1 ]; then
 if [ ! -d $procdir/geo_500m ]; then
  echo "generating low resolution geo files (500 m)"
  thisdir=`pwd`
  cd $procdir
  submit_geo_lores.py
  geodir=geo_500m
 fi
 if [ -f $procdir/geo_500m/locked ]; then
  echo "seems like the lores geocoding is locked by another process. please wait for it to finish first"
  exit
 fi

 GEOCDIR=GEOC_500m
 mkdir -p $GEOCDIR
 cd $thisdir
fi

master=`ls $procdir/$geodir/*[0-9].hgt | rev | cut -d '/' -f1 | rev | cut -d '.' -f1 | head -n1`

echo "Processing dir: $procdir"
echo "Master image: $master"
echo "Interferogram: $ifg"

width=`awk '$1 == "range_samples:" {print $2}' ${procdir}/RSLC/$master/$master.rslc.mli.par`; #echo $width
length=`awk '$1 == "azimuth_lines:" {print $2}' ${procdir}/RSLC/$master/$master.rslc.mli.par`; #echo $length
reducfac=`echo $width | awk '{if(int($1/16000) > 1) print int($1/16000); else print 1}'`

#if [ -d "${procdir}/GEOC" ]; then rm -rf ${procdir}/GEOC; fi
if [ ! -d "${procdir}/$GEOCDIR" ]; then mkdir ${procdir}/$GEOCDIR 2>/dev/null; fi

lat=`awk '$1 == "corner_lat:" {print $2}' ${procdir}/$geodir/EQA.dem_par`
lon=`awk '$1 == "corner_lon:" {print $2}' ${procdir}/$geodir/EQA.dem_par`
latstep=`awk '$1 == "post_lat:" {print $2}' ${procdir}/$geodir/EQA.dem_par`
lonstep=`awk '$1 == "post_lon:" {print $2}' ${procdir}/$geodir/EQA.dem_par`
length_dem=`awk '$1 == "nlines:" {print $2}' ${procdir}/$geodir/EQA.dem_par`
width_dem=`awk '$1 == "width:" {print $2}' ${procdir}/$geodir/EQA.dem_par`
reducfac_dem=`echo $width_dem | awk '{if(int($1/2000) > 1) print int($1/2000); else print 1}'`
lat1=`echo $lat $latstep $length_dem | awk '{printf "%7f", ($1+($2*($3-1)));}'` # Substract one because width starts at zero
lon1=`echo $lon $lonstep $width_dem | awk '{printf "%7f", ($1+($2*($3-1)));}'`
latstep=`awk '$1 == "post_lat:" {if($2<0) printf "%7f", -1*$2; else printf "%7f", $2}' ${procdir}/$geodir/EQA.dem_par`
lonstep=`awk '$1 == "post_lon:" {if($2<0) printf "%7f", -1*$2; else printf "%7f", $2}' ${procdir}/$geodir/EQA.dem_par`
# Because no wavelength is reported in master.rmli.par file, we calculated here according to the radar frequency (IN CENTIMETERS)
# Frequency = (C / Wavelength), Where: Frequency: Frequency of the wave in hertz (hz). C: Speed of light (29,979,245,800 cm/sec (3 x 10^10 approx))
lambda=`awk '$1 == "radar_frequency:" {print 29979245800/$2}' ${procdir}/RSLC/$master/$master.rslc.mli.par`;

echo " Running doGeocoding step for unwrapped"
#this is for the case if we start the geotiff generation within or outside of licsar_make_frame:
if [ -d ${procdir}/LOGS ]; then
 logfile=${procdir}/LOGS/13_doGeocoding_$ifg.log
else
 logfile=${procdir}/13_doGeocoding_$ifg.log
fi
echo "   check "$logfile" if something goes wrong "
#rm -f $logfile

#mdate=`echo $ifg | awk '{print $2}'`;
#sdate=`echo $ifg | awk '{print $3}'`;
#if [ -d "${procdir}/GEOC/${ifg}" ]; then rm -rf ${procdir}/GEOC/${ifg}; fi
if [ ! -d "${procdir}/$GEOCDIR/${ifg}" ]; then mkdir ${procdir}/$GEOCDIR/${ifg}; fi
echo "   Geocoding results for inteferogram: ${ifg}" ;

if [ $UNWDO -eq 1 ]; then
# Unwrapped interferogram
if [ -e ${procdir}/IFG/${ifg}/${ifg}.unw ]; then
  # Geocode all the data
  if [ ! -e ${procdir}/$GEOCDIR/${ifg}/${ifg}.geo.unw ]; then
  echo "Replacing nan for 0 values (trick by Yu Morishita)"
   #I will correct nan-corrected unwrapped ifg here and then convert the unwrapped phase to metric units
   replace_values ${procdir}/IFG/${ifg}/${ifg}.unw nan 0 ${procdir}/IFG/${ifg}/${ifg}.unw0 $width >> $logfile
   if [ `ls -al ${procdir}/IFG/${ifg}/${ifg}.unw  | gawk {'print $5'}` -eq  `ls -al ${procdir}/IFG/${ifg}/${ifg}.unw0  | gawk {'print $5'}` ]; then
     mv ${procdir}/IFG/${ifg}/${ifg}.unw0 ${procdir}/IFG/${ifg}/${ifg}.unw
   else
     rm ${procdir}/IFG/${ifg}/${ifg}.unw0 
   fi
   #mv ${procdir}/IFG/${ifg}/${ifg}.disp ${procdir}/IFG/${ifg}/${ifg}.unw
   echo "Geocoding"
   geocode_back ${procdir}/IFG/${ifg}/${ifg}.unw $width ${procdir}/$geodir/$master.lt_fine ${procdir}/$GEOCDIR/${ifg}/${ifg}.geo.unw ${width_dem} ${length_dem} 0 0 >> $logfile
  fi
#  if [ ! -e ${procdir}/GEOC/${ifg}/${ifg}.geo.disp ]; then
#   echo "Converting unwrapped phase to displacements"
#   dispmap ${procdir}/IFG/${ifg}/${ifg}.unw ${procdir}/$geodir/$master.hgt ${procdir}/SLC/$master/$master.slc.par ${procdir}/IFG/${ifg}/${ifg}.off ${procdir}/IFG/${ifg}/${ifg}.disp 0 0 >> $logfile
#   geocode_back ${procdir}/IFG/${ifg}/${ifg}.disp $width ${procdir}/$geodir/$master.lt_fine ${procdir}/GEOC/${ifg}/${ifg}.geo.disp ${width_dem} ${length_dem} 0 0 >> $logfile   
#  fi
  # Convert to geotiff
  if [ ! -e ${procdir}/$GEOCDIR/${ifg}/${ifg}.geo.unw.tif ]; then 
   echo "Converting to GeoTIFF"
   data2geotiff ${procdir}/$geodir/EQA.dem_par ${procdir}/$GEOCDIR/${ifg}/${ifg}.geo.unw 2 ${procdir}/$GEOCDIR/${ifg}/${ifg}.geo.unw.orig.tif 0.0  >> $logfile 2>/dev/null
   gdal_translate -of GTiff -ot Float32 -co COMPRESS=DEFLATE -co PREDICTOR=3 ${procdir}/$GEOCDIR/${ifg}/${ifg}.geo.unw.orig.tif ${procdir}/$GEOCDIR/${ifg}/${ifg}.geo.unw.tif >> $logfile 2>/dev/null
   rm ${procdir}/$GEOCDIR/${ifg}/${ifg}.geo.unw.orig.tif
#   data2geotiff ${procdir}/$geodir/EQA.dem_par ${procdir}/GEOC/${ifg}/${ifg}.geo.disp 2 ${procdir}/GEOC/${ifg}/${ifg}.geo.disp.tif 0.0  >> $logfile 2>/dev/null
  fi
  # Create bmps
  #ras_linear ${procdir}/GEOC/${ifg}/${ifg}.geo.unw ${width_dem} - - $reducfac_dem $reducfac_dem - - - ${procdir}/GEOC/${ifg}/${ifg}.geo.unw.bmp >> $logfile
  if [ ! -e ${procdir}/$GEOCDIR/${ifg}/${ifg}.geo.unw_blk.bmp ]; then
   echo "Converting to raster previews"
   rasrmg ${procdir}/$GEOCDIR/${ifg}/${ifg}.geo.unw ${procdir}/$geodir/EQA.${master}.slc.mli ${width_dem} - - - $reducfac_dem $reducfac_dem - - - - - ${procdir}/$GEOCDIR/${ifg}/${ifg}.geo.unw_blk.bmp >> $logfile
  fi 
   #rasrmg ${procdir}/GEOC/${ifg}/${ifg}.geo.disp ${procdir}/$geodir/EQA.${master}.slc.mli ${width_dem} - - - $reducfac_dem $reducfac_dem - - - - - ${procdir}/GEOC/${ifg}/${ifg}.geo.disp_blk.bmp >> $logfile
   #get min and max for disp image
#   gdalinfo -stats ${procdir}/GEOC/${ifg}/${ifg}.geo.disp.tif > ${procdir}/GEOC/${ifg}/${ifg}.geo.disp.stats
#   min=`grep MINIMUM ${procdir}/GEOC/${ifg}/${ifg}.geo.disp.stats | cut -d '=' -f2`
#   max=`grep MAXIMUM ${procdir}/GEOC/${ifg}/${ifg}.geo.disp.stats | cut -d '=' -f2`
#   std=`grep STDDEV ${procdir}/GEOC/${ifg}/${ifg}.geo.disp.stats | cut -d '=' -f2`
#   max=`echo $max-$std | bc`
#   min=`echo $min+$std | bc`
   
   #generate displacement image
#   visdt_pwr.py ${procdir}/GEOC/${ifg}/${ifg}.geo.disp ${procdir}/$geodir/EQA.${master}.slc.mli ${width_dem} $min $max -b -p ${procdir}/GEOC/${ifg}/${ifg}.geo.disp_blk.png 2>/dev/null
#   convert ${procdir}/GEOC/${ifg}/${ifg}.geo.disp_blk.png ${procdir}/GEOC/${ifg}/${ifg}.geo.disp_blk.bmp
   #rasdt_cmap ${procdir}/GEOC/${ifg}/${ifg}.geo.disp ${procdir}/$geodir/EQA.${master}.slc.mli ${width_dem} - - - $reducfac_dem $reducfac_dem $min $max 0 - - - ${procdir}/GEOC/${ifg}/${ifg}.geo.disp_blk.bmp >> $logfile   
  if [ ! -e ${procdir}/$GEOCDIR/${ifg}/${ifg}.geo.unw.bmp ]; then
   convert -transparent black -resize $RESIZE'%' ${procdir}/$GEOCDIR/${ifg}/${ifg}.geo.unw_blk.bmp ${procdir}/$GEOCDIR/${ifg}/${ifg}.geo.unw.png
  fi
#  convert -transparent black ${procdir}/GEOC/${ifg}/${ifg}.geo.disp_blk.bmp ${procdir}/GEOC/${ifg}/${ifg}.geo.disp.png
fi
fi

#if [ $UNFILT -eq 1 ]; then
# ifgext="diff"
# ifgtext="unfiltered"
#else
# ifgext="filt.diff"
# ifgtext="filtered"
#fi
#echo "bagr"

if [ $IFGDO -eq 1 ]; then
ifgext="filt.diff"
ifgtext="filtered"
ifgout="diff"
if [ -e ${procdir}/IFG/${ifg}/${ifg}.$ifgext ] && [ ! -e ${procdir}/$GEOCDIR/${ifg}/${ifg}.geo.diff_pha.tif ]; then
 echo "Creating "$ifgtext" interferogram tiffs"
 # Extract the mag and phase
 #cpx_to_real ${procdir}/IFG/${ifg}/${ifg}.filt.diff ${procdir}/IFG/${ifg}/${ifg}.diff_mag $width 3  >> $logfile
 cpx_to_real ${procdir}/IFG/${ifg}/${ifg}.$ifgext ${procdir}/IFG/${ifg}/${ifg}.$ifgout'_pha' $width 4  >> $logfile
 # Geocode all the data
 geocode_back ${procdir}/IFG/${ifg}/${ifg}.$ifgext $width ${procdir}/$geodir/$master.lt_fine ${procdir}/$GEOCDIR/${ifg}/${ifg}.geo.$ifgout ${width_dem} ${length_dem} 1 1 >> $logfile
 #geocode_back ${procdir}/IFG/${ifg}/${ifg}.diff_mag $width ${procdir}/$geodir/$master.lt_fine ${procdir}/GEOC/${ifg}/${ifg}.geo.diff_mag ${width_dem} ${length_dem} 1 0 >> $logfile
 geocode_back ${procdir}/IFG/${ifg}/${ifg}.$ifgout'_pha' $width ${procdir}/$geodir/$master.lt_fine ${procdir}/$GEOCDIR/${ifg}/${ifg}.geo.$ifgout'_pha' ${width_dem} ${length_dem} 0 0 >> $logfile
 # Convert to geotiff
 #data2geotiff ${procdir}/$geodir/EQA.dem_par ${procdir}/GEOC/${ifg}/${ifg}.geo.diff_mag 2 ${procdir}/GEOC/${ifg}/${ifg}.geo.diff_mag.tif 0.0  >> $logfile 2>/dev/null
 data2geotiff ${procdir}/$geodir/EQA.dem_par ${procdir}/$GEOCDIR/${ifg}/${ifg}.geo.$ifgout'_pha' 2 ${procdir}/$GEOCDIR/${ifg}/${ifg}.geo.$ifgout'_pha.orig.tif' 0.0  >> $logfile 2>/dev/null
 gdal_translate -of GTiff -ot Float32 -co COMPRESS=DEFLATE -co PREDICTOR=3 ${procdir}/$GEOCDIR/${ifg}/${ifg}.geo.$ifgout'_pha.orig.tif' ${procdir}/$GEOCDIR/${ifg}/${ifg}.geo.$ifgout'_pha.tif' >> $logfile 2>/dev/null
 rm ${procdir}/$GEOCDIR/${ifg}/${ifg}.geo.$ifgout'_pha.orig.tif'
 # Create bmps
 rasmph_pwr ${procdir}/$GEOCDIR/${ifg}/${ifg}.geo.$ifgout ${procdir}/$geodir/EQA.${master}.slc.mli ${width_dem} - - - $reducfac_dem $reducfac_dem - - - ${procdir}/$GEOCDIR/${ifg}/${ifg}.geo.$ifgout'_blk.bmp' >> $logfile
 convert -transparent black -resize $RESIZE'%' ${procdir}/$GEOCDIR/${ifg}/${ifg}.geo.$ifgout'_blk.bmp' ${procdir}/$GEOCDIR/${ifg}/${ifg}.geo.$ifgout.png
 #rasmph_pwr ${procdir}/GEOC/${ifg}/${ifg}.geo.diff_mag ${procdir}/$geodir/EQA.${master}.slc.mli ${width_dem} - - - $reducfac_dem $reducfac_dem - - - ${procdir}/GEOC/${ifg}/${ifg}.geo.diff_mag.bmp >> $logfile
 #rasmph_pwr ${procdir}/GEOC/${ifg}/${ifg}.geo.diff_pha ${procdir}/$geodir/EQA.${master}.slc.mli ${width_dem} - - - $reducfac_dem $reducfac_dem - - - ${procdir}/GEOC/${ifg}/${ifg}.geo.diff_pha.bmp >> $logfile

 #This was to generate only mag and pha bmp..
 #raspwr ${procdir}/GEOC/${ifg}/${ifg}.geo.diff_mag ${width_dem} - - $reducfac_dem $reducfac_dem - - - ${procdir}/GEOC/${ifg}/${ifg}.geo.diff_mag.bmp >> $logfile
 #ras_linear ${procdir}/GEOC/${ifg}/${ifg}.geo.diff_pha ${width_dem} - - $reducfac_dem $reducfac_dem - - - ${procdir}/GEOC/${ifg}/${ifg}.geo.diff_pha.bmp >> $logfile
 # Clean
 rm ${procdir}/IFG/${ifg}/${ifg}.$ifgout'_mag' ${procdir}/IFG/${ifg}/${ifg}.$ifgout'_pha' 2>/dev/null
fi
fi

if [ $UNFILT -eq 1 ]; then
 #do the unfiltered ifg as well as filtered
 ifgtext="unfiltered"
 ifgext="diff"
 ifgout="diff_unfiltered"
 echo "Creating "$ifgtext" interferogram tiffs"
 cpx_to_real ${procdir}/IFG/${ifg}/${ifg}.$ifgext ${procdir}/IFG/${ifg}/${ifg}.$ifgout'_pha' $width 4  >> $logfile
 geocode_back ${procdir}/IFG/${ifg}/${ifg}.$ifgext $width ${procdir}/$geodir/$master.lt_fine ${procdir}/$GEOCDIR/${ifg}/${ifg}.geo.$ifgout ${width_dem} ${length_dem} 1 1 >> $logfile
 geocode_back ${procdir}/IFG/${ifg}/${ifg}.$ifgout'_pha' $width ${procdir}/$geodir/$master.lt_fine ${procdir}/$GEOCDIR/${ifg}/${ifg}.geo.$ifgout'_pha' ${width_dem} ${length_dem} 0 0 >> $logfile
 data2geotiff ${procdir}/$geodir/EQA.dem_par ${procdir}/$GEOCDIR/${ifg}/${ifg}.geo.$ifgout'_pha' 2 ${procdir}/$GEOCDIR/${ifg}/${ifg}.geo.$ifgout'_pha.orig.tif' 0.0  >> $logfile 2>/dev/null
 gdal_translate -of GTiff -ot Float32 -co COMPRESS=DEFLATE -co PREDICTOR=3 ${procdir}/$GEOCDIR/${ifg}/${ifg}.geo.$ifgout'_pha.orig.tif' ${procdir}/$GEOCDIR/${ifg}/${ifg}.geo.$ifgout'_pha.tif' >> $logfile 2>/dev/null
 rm ${procdir}/$GEOCDIR/${ifg}/${ifg}.geo.$ifgout'_pha.orig.tif'
 rasmph_pwr ${procdir}/$GEOCDIR/${ifg}/${ifg}.geo.$ifgout ${procdir}/$geodir/EQA.${master}.slc.mli ${width_dem} - - - $reducfac_dem $reducfac_dem - - - ${procdir}/$GEOCDIR/${ifg}/${ifg}.geo.$ifgout'_blk.bmp' >> $logfile
 convert -transparent black -resize $RESIZE'%' ${procdir}/$GEOCDIR/${ifg}/${ifg}.geo.$ifgout'_blk.bmp' ${procdir}/$GEOCDIR/${ifg}/${ifg}.geo.$ifgout.png
 rm ${procdir}/IFG/${ifg}/${ifg}.$ifgout'_mag' ${procdir}/IFG/${ifg}/${ifg}.$ifgout'_pha' 2>/dev/null
fi


if [ $COHDO -eq 1 ]; then
# Unfiltered coherence
if [ -e ${procdir}/IFG/${ifg}/${ifg}.cc ] && [ ! -e ${procdir}/$GEOCDIR/${ifg}/${ifg}.geo.cc.tif ]; then
  echo "Creating unfiltered coherence tiffs"
  # Geocode
  geocode_back ${procdir}/IFG/${ifg}/${ifg}.cc $width ${procdir}/$geodir/$master.lt_fine ${procdir}/$GEOCDIR/${ifg}/${ifg}.geo.cc ${width_dem} ${length_dem} 1 0 >> $logfile
  # Convert to geotiff
  data2geotiff ${procdir}/$geodir/EQA.dem_par ${procdir}/$GEOCDIR/${ifg}/${ifg}.geo.cc 2 ${procdir}/$GEOCDIR/${ifg}/${ifg}.geo.cc.orig.tif 0.0  >> $logfile 2>/dev/null
  #for compression types differences, check e.g. https://kokoalberti.com/articles/geotiff-compression-optimization-guide/
  gdal_translate -of GTiff -ot Byte -scale 0 1 0 255 -co COMPRESS=DEFLATE -co PREDICTOR=2 ${procdir}/$GEOCDIR/${ifg}/${ifg}.geo.cc.orig.tif ${procdir}/$GEOCDIR/${ifg}/${ifg}.geo.cc.tif >> $logfile 2>/dev/null
  rm ${procdir}/$GEOCDIR/${ifg}/${ifg}.geo.cc.orig.tif
  # create bmps
  #new version of gamma shows coherence in colour... using old-school cpxfiddle as workaround
  #rascc ${procdir}/GEOC/${ifg}/${ifg}.geo.cc - ${width_dem} - - - $reducfac_dem $reducfac_dem 0 1 - - - ${procdir}/GEOC/${ifg}/${ifg}.geo.cc_blk.bmp >> $logfile
 # module load doris
  byteswap -o ${ifg}.geo.cc.tmp 4 ${procdir}/$GEOCDIR/${ifg}/${ifg}.geo.cc >/dev/null
  cpxfiddle -q normal -w ${width_dem} -f r4 -o sunraster -c gray -M $reducfac_dem/$reducfac_dem -r 0/0.9 ${ifg}.geo.cc.tmp > ${procdir}/$GEOCDIR/${ifg}/${ifg}.geo.cc_blk.ras 2>/dev/null
  convert ${procdir}/$GEOCDIR/${ifg}/${ifg}.geo.cc_blk.ras ${procdir}/$GEOCDIR/${ifg}/${ifg}.geo.cc_blk.bmp
  rm -f ${procdir}/$GEOCDIR/${ifg}/${ifg}.geo.cc_blk.ras ${ifg}.geo.cc.tmp
  # Need to remove the black border, but the command below is no good as it removes the black parts of the coherence!
  #convert -transparent black -resize 30% ${procdir}/GEOC/${ifg}/${ifg}.geo.cc_blk.bmp ${procdir}/GEOC/${ifg}/${ifg}.geo.cc.bmp
  convert -transparent black -resize $RESIZE'%' ${procdir}/$GEOCDIR/${ifg}/${ifg}.geo.cc_blk.bmp ${procdir}/$GEOCDIR/${ifg}/${ifg}.geo.cc.png
fi
fi

if [ $MLIDO -eq 1 ]; then
#finally amplitudes/intensities
mkdir $procdir/GEOC.MLI 2>/dev/null
for im in `echo $ifg | sed 's/_/ /'`; do
if [ -e ${procdir}/RSLC/$im/$im.rslc.mli ] && [ ! -d ${procdir}/GEOC.MLI/$im ]; then
 mkdir -p ${procdir}/GEOC.MLI/$im
 echo "preparing GeoTIFF and PNG for intensity image of "$im
 #geocode MLI
 geocode_back ${procdir}/RSLC/$im/$im.rslc.mli $width ${procdir}/$geodir/$master.lt_fine ${procdir}/GEOC.MLI/$im/$im.geo.mli ${width_dem} ${length_dem} 0 0 >> $logfile
 #convert MLI to geotiff
 data2geotiff ${procdir}/$geodir/EQA.dem_par ${procdir}/GEOC.MLI/$im/$im.geo.mli 2 ${procdir}/GEOC.MLI/$im/$im.geo.mli.tif 0.0  >> $logfile 2>/dev/null
 #generate raster preview
 raspwr ${procdir}/GEOC.MLI/$im/$im.geo.mli ${width_dem} - - $reducfac_dem $reducfac_dem - - - ${procdir}/GEOC.MLI/$im/$im.geo.mli.bmp 0 - >> $logfile
 convert -transparent black -resize $RESIZE'%' ${procdir}/GEOC.MLI/$im/$im.geo.mli.bmp ${procdir}/GEOC.MLI/$im/$im.geo.mli.png
fi 
done
fi

echo "done"

#Move it to the public area(?)
#for filename in .......; do
# if [ ! -f ${publicdir}/ ]
