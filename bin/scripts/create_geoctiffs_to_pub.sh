#!/bin/bash

source $LiCSARpath/lib/LiCSAR_bash_lib.sh

if [ -z $2 ]; then 
 echo "inputs are: create_geoctiffs_to_pub.sh procdir ifg, e.g. \`pwd\` 20160101_20160202";
 echo "optional parameters:"
 #echo "-c lat1/lat2/lon1/lon2 - would establish a crop area for the frame"
 echo "-u .... geocode also unfiltered wrapped interferogram"
 echo "-F .... do full resolution previews" # (needed for KML) - this will use different colour bar"
 echo "-L .... do low resolution geotiffs (500x500 m)"
 echo "-H .... do full resolution geotiffs (50x50 m)"
 echo "-m .... use land mask if available (we turn this on by default now)"
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
FULL=0
mask=1
clip=0
magcc=1

while getopts ":uabmUCFHIML" option; do
 case "${option}" in
  a) MLIDO=0; echo "auto-regime - doing all except of MLIs";
     ;;
  b) UNFILT=0; IFGDO=0; COHDO=0; UNWDO=0; MLIDO=0; echo "bypassing other (fast fix)";
     ;;
  u) UNFILT=1; echo "unfiltered ifg will be geocoded";
     ;;
  m) mask=1; echo "will try masking"
     ;;
  U) IFGDO=0; COHDO=0; MLIDO=0; echo "you do ONLY unw geo now..";
     ;;
  C) IFGDO=0; UNWDO=0; MLIDO=0; echo "you do ONLY coh geo now..";
     ;;
  I) COHDO=0; UNWDO=0; MLIDO=0; IFGDO=1; echo "you do ONLY ifg geo now..";
     ;;
  M) COHDO=0; UNWDO=0; MLIDO=1; IFGDO=0; echo "you do ONLY MLIs geo now..";
     ;;
  F) FULL=1; echo "previews will be generated in full resolution";
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
 fi
 if [ -f $procdir/geo_50m/locked ]; then
  echo "seems like the hires geocoding is locked by another process. please wait for it to finish first"
  exit
 fi
 geodir=geo_50m
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
 fi
 if [ -f $procdir/geo_500m/locked ]; then
  echo "seems like the lores geocoding is locked by another process. please wait for it to finish first"
  exit
 fi
 geodir=geo_500m
 GEOCDIR=GEOC_500m
 mkdir -p $GEOCDIR
 cd $thisdir
fi

master=`ls $procdir/$geodir/*[0-9].hgt | rev | cut -d '/' -f1 | rev | cut -d '.' -f1 | head -n1`
#make sure we use abs paths...
procdir=`realpath $procdir`

#this is just to check if master is properly in RSLCs...:
if [ ! -e ${procdir}/RSLC/$master/$master.rslc.mli ]; then 
 if [ -e ${procdir}/SLC/$master/$master.slc.mli ]; then 
  echo "correcting bad links for master file"; 
  rm -r ${procdir}/RSLC/$master
  mkdir ${procdir}/RSLC/$master
  cd ${procdir}/RSLC/$master
  for x in `ls ${procdir}/SLC/$master/* -d`; do
    ln -s $x `basename $x | sed 's/\.slc\./\.rslc\./'`
  done
  chmod -R 777 ${procdir}/RSLC/$master
 fi
fi


echo "Processing dir: $procdir"
echo "Master image: $master"
echo "Interferogram: $ifg"

width=`awk '$1 == "range_samples:" {print $2}' ${procdir}/SLC/$master/$master.slc.mli.par`; #echo $width
length=`awk '$1 == "azimuth_lines:" {print $2}' ${procdir}/SLC/$master/$master.slc.mli.par`; #echo $length
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
lambda=`awk '$1 == "radar_frequency:" {print 29979245800/$2}' ${procdir}/SLC/$master/$master.slc.mli.par`;

# update 2022: set hgt file and warp towards it
frame=`basename $procdir`
tr=`track_from_frame $frame`
hgtfile=$LiCSAR_public/$tr/$frame/metadata/$frame.geo.hgt.tif
if [ -z $hgtfile ]; then
 hgtfile=`ls GEOC/geo/*.geo.hgt.tif 2>/dev/null | head -n 1 2>/dev/null`
fi

#echo "Running doGeocoding step" #" for unwrapped"
#this is for the case if we start the geotiff generation within or outside of licsar_make_frame:
if [ -d ${procdir}/LOGS ]; then
 logfile=${procdir}/LOGS/13_doGeocoding_$ifg.log
else
 logfile=${procdir}/13_doGeocoding_$ifg.log
fi
#echo "check "$logfile" if something goes wrong "
#rm -f $logfile

#mdate=`echo $ifg | awk '{print $2}'`;
#sdate=`echo $ifg | awk '{print $3}'`;
#if [ -d "${procdir}/GEOC/${ifg}" ]; then rm -rf ${procdir}/GEOC/${ifg}; fi

#echo "   Geocoding results for inteferogram: ${ifg}" ;

if [ $UNWDO -eq 1 ]; then
if [ ! -d "${procdir}/$GEOCDIR/${ifg}" ]; then mkdir ${procdir}/$GEOCDIR/${ifg}; fi
# Unwrapped interferogram
if [ -f ${procdir}/IFG/${ifg}/${ifg}.unw ]; then
  # Geocode all the data
if [ ! -e ${procdir}/$GEOCDIR/${ifg}/${ifg}.geo.unw ]; then
  echo "geocoding unwrapped data"
  #echo "Replacing nan for 0 values (trick by Yu Morishita)"
   #I will correct nan-corrected unwrapped ifg here and then convert the unwrapped phase to metric units
   replace_values ${procdir}/IFG/${ifg}/${ifg}.unw nan 0 ${procdir}/IFG/${ifg}/${ifg}.unw0 $width >> $logfile
   if [ `ls -al ${procdir}/IFG/${ifg}/${ifg}.unw  | gawk {'print $5'}` -eq  `ls -al ${procdir}/IFG/${ifg}/${ifg}.unw0  | gawk {'print $5'}` ]; then
     mv ${procdir}/IFG/${ifg}/${ifg}.unw0 ${procdir}/IFG/${ifg}/${ifg}.unw
   else
     rm ${procdir}/IFG/${ifg}/${ifg}.unw0 
   fi
   #mv ${procdir}/IFG/${ifg}/${ifg}.disp ${procdir}/IFG/${ifg}/${ifg}.unw
   #echo "Geocoding"
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
   #shifting by median
   nctempfile=${procdir}/$GEOCDIR/${ifg}/${ifg}.geo.unw.orig.nc
   gmt grdclip ${procdir}/$GEOCDIR/${ifg}/${ifg}.geo.unw.orig.tif -G$nctempfile -Sr0/NaN
   gmt grdmath $nctempfile $nctempfile MEDIAN SUB = ${procdir}/$GEOCDIR/${ifg}/${ifg}.geo.unw.orig.tif=gd:gtiff
   minmaxcolour=`gmt grdinfo -T+a0.1+s $nctempfile`
   minmaxreal=`gmt grdinfo -T $nctempfile`
   # update 2022: align to hgt - should work also for clipping!
   if [ -f $hgtfile ]; then
    gdalwarp2match.py ${procdir}/$GEOCDIR/${ifg}/${ifg}.geo.unw.orig.tif $hgtfile ${procdir}/$GEOCDIR/${ifg}/${ifg}.geo.unw.orig2.tif
    mv ${procdir}/$GEOCDIR/${ifg}/${ifg}.geo.unw.orig2.tif ${procdir}/$GEOCDIR/${ifg}/${ifg}.geo.unw.orig.tif
   fi
   gdal_translate -of GTiff -ot Float32 -co COMPRESS=DEFLATE -co PREDICTOR=3 -a_srs EPSG:4326 ${procdir}/$GEOCDIR/${ifg}/${ifg}.geo.unw.orig.tif ${procdir}/$GEOCDIR/${ifg}/${ifg}.geo.unw.tif >> $logfile 2>/dev/null
   rm ${procdir}/$GEOCDIR/${ifg}/${ifg}.geo.unw.orig.tif
#   data2geotiff ${procdir}/$geodir/EQA.dem_par ${procdir}/GEOC/${ifg}/${ifg}.geo.disp 2 ${procdir}/GEOC/${ifg}/${ifg}.geo.disp.tif 0.0  >> $logfile 2>/dev/null
   echo "Generating preview PNG"
  if [ $mask -eq 1 ]; then
    echo "testing - include landmask (should work if procdir ends by frame name)"
    #exit
    frame=`basename $procdir`
    # source $LiCSARpath/lib/LiCSAR_bash_lib.sh
    create_preview_unwrapped ${procdir}/$GEOCDIR/${ifg}/${ifg}.geo.unw.tif $frame
  else
       unw_bmp=${procdir}/$GEOCDIR/${ifg}/${ifg}.geo.unw.png
       scalebar_bmp=${procdir}/$GEOCDIR/${ifg}/${ifg}.geo.unw.scale.png
       unwcpt=${procdir}/$GEOCDIR/${ifg}/unw.cpt
       gmt makecpt -C$LiCSARpath/misc/colourmap.cpt -Iz $minmaxcolour/0.025 >$unwcpt
       #frame=`echo ${procdir} | rev | cut -d '/' -f1 | rev`
       #tr=`echo $frame | cut -c -3 | sed 's/^0//' | sed 's/^0//'`
       #hgtfile=$LiCSAR_public/$tr/$frame/metadata/$frame.geo.hgt.tif
       hillshadecmd=''
       hillshadefile=''
       if [ -f $hgtfile ]; then
         hillshadefile=$LiCSAR_public/$tr/$frame/metadata/$frame.geo.hillshade.nc
        if [ ! -f $hillshadefile ]; then
          hillshadefile=${procdir}/$GEOCDIR/${ifg}/hillshade.nc
          echo "generating hillshade"
        opass=`echo $frame | cut -c 4`
        if [ $opass == 'A' ]; then
          deg=258
        else
          deg=102
        fi
        gmt grdgradient -A$deg -Nt1 -G$hillshadefile $hgtfile
        fi
       else
         echo "hgt geofile not found. will skip generating hillshade"
       fi
       if [ -f $hillshadefile ]; then
        hillshadecmd="-I"$hillshadefile
       fi
      
       gmt grdimage $nctempfile -C$unwcpt -JM1 -Q -nn+t0.1 $hillshadecmd -A$unw_bmp 2>/dev/null
       if [ ! -f $unw_bmp ]; then
        echo "some error, perhaps with hillshade.. skipping it"
        gmt grdimage $nctempfile -C$unwcpt -JM1 -Q -nn+t0.1 -A$unw_bmp.tt.png
        #in such case we should be ok with 255 colour palette
        convert $unw_bmp.tt.png PNG8:$unw_bmp
        rm $unw_bmp.tt.png
       fi

      #need to prepare a colorbar based on these values!!!!
      #you know... the rounding here is not really important... just a colourbar.. or not?
      mincol=`echo $minmaxcolour | cut -d '/' -f1 | cut -d 'T' -f2`
      maxcol=`echo $minmaxcolour | cut -d '/' -f2 `
      #expecting sentinel
      minval=`python -c 'print(round('$mincol'*5.546/(4*3.14159265)))'`
      maxval=`python -c 'print(round('$maxcol'*5.546/(4*3.14159265)))'`
      #add also real min and max values
      minreal=`echo $minmaxreal | cut -d 'T' -f2 | cut -d '/' -f1 | cut -d '.' -f1`
      maxreal=`echo $minmaxreal | cut -d '/' -f2 | cut -d '.' -f1`
      minrealval=`python -c 'print(round('$minreal'*5.546/(4*3.14159265)))'`
      maxrealval=`python -c 'print(round('$maxreal'*5.546/(4*3.14159265)))'`
      #burn them to the scalebar
      minvalsize=`echo $minval | wc -m `
      if [ $minvalsize -gt 4 ]; then
       xsize=20
      elif [ $minvalsize -eq 4 ]; then
       xsize=40
      elif [ $minvalsize -eq 3 ]; then
       xsize=60
      else
       xsize=80
      fi
      convert -font helvetica -fill black -pointsize 40 -draw "text "$xsize",115 '"$minval"'" $LiCSARpath/misc/scalebar_unwrapped_empty.png $scalebar_bmp.temp.png
      convert -font helvetica -fill black -pointsize 40 -draw "text 1100,115 '"$maxval" cm'" $scalebar_bmp.temp.png $scalebar_bmp
      mv $scalebar_bmp $scalebar_bmp.temp.png
      #add real values
      convert -font helvetica -fill black -pointsize 35 -draw "text "$xsize",165 '[min "$minrealval" cm]'" $scalebar_bmp.temp.png  $scalebar_bmp
      mv $scalebar_bmp $scalebar_bmp.temp.png
      convert -font helvetica -fill black -pointsize 35 -draw "text 1020,165 '[max "$maxrealval" cm]'" $scalebar_bmp.temp.png  $scalebar_bmp

      convert $unw_bmp -resize 680x \( $scalebar_bmp -resize 400x  -background none -gravity center \) -gravity southwest -geometry +7+7 -composite -flatten -transparent black $unw_bmp.sm.png
      if [ $FULL -eq 1 ]; then
       mv $unw_bmp ${procdir}/$GEOCDIR/${ifg}/${ifg}.geo.unw.full.png
      fi
      mv $unw_bmp.sm.png $unw_bmp
      
      rm $unwcpt $nctempfile $scalebar_bmp.temp.png ${procdir}/$GEOCDIR/${ifg}/hillshade.nc 2>/dev/null
      rm ${procdir}/$GEOCDIR/${ifg}/${ifg}.geo.unw.orig.tif $scalebar_bmp 2>/dev/null
  fi
fi


fi










  #ras_linear ${procdir}/GEOC/${ifg}/${ifg}.geo.unw ${width_dem} - - $reducfac_dem $reducfac_dem - - - ${procdir}/GEOC/${ifg}/${ifg}.geo.unw.bmp >> $logfile
  #if [ ! -e ${procdir}/$GEOCDIR/${ifg}/${ifg}.geo.unw_blk.bmp ]; then
  # echo "Converting to raster previews"
  # rasrmg ${procdir}/$GEOCDIR/${ifg}/${ifg}.geo.unw ${procdir}/$geodir/EQA.${master}.slc.mli ${width_dem} - - - $reducfac_dem $reducfac_dem - - - - - ${procdir}/$GEOCDIR/${ifg}/${ifg}.geo.unw_blk.bmp >> $logfile
  #fi 
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
  #if [ ! -e ${procdir}/$GEOCDIR/${ifg}/${ifg}.geo.unw.bmp ]; then
  # convert -transparent black -resize $RESIZE'%' ${procdir}/$GEOCDIR/${ifg}/${ifg}.geo.unw_blk.bmp ${procdir}/$GEOCDIR/${ifg}/${ifg}.geo.unw.png
   #if [ $FULL -eq 1 ]; then
   # convert -transparent black ${procdir}/$GEOCDIR/${ifg}/${ifg}.geo.unw_blk.bmp ${procdir}/$GEOCDIR/${ifg}/${ifg}.geo.unw.full.png
   #fi
  #fi
  #rm ${procdir}/$GEOCDIR/${ifg}/${ifg}.geo.unw_blk.bmp 2>/dev/null
#  convert -transparent black ${procdir}/GEOC/${ifg}/${ifg}.geo.disp_blk.bmp ${procdir}/GEOC/${ifg}/${ifg}.geo.disp.png
#fi
fi  #end of unwdo



#if [ $UNFILT -eq 1 ]; then
# ifgext="diff"
# ifgtext="unfiltered"
#else
# ifgext="filt.diff"
# ifgtext="filtered"
#fi
#echo "bagr"

if [ $IFGDO -eq 1 ]; then
if [ ! -d "${procdir}/$GEOCDIR/${ifg}" ]; then mkdir ${procdir}/$GEOCDIR/${ifg}; fi
ifgext="filt.diff"
ifgtext="filtered"
ifgout="diff"
if [ -e ${procdir}/IFG/${ifg}/${ifg}.$ifgext ] && [ ! -e ${procdir}/$GEOCDIR/${ifg}/${ifg}.geo.diff_pha.tif ]; then
 echo "Creating "$ifgtext" interferogram tiffs"
 # Extract the mag and phase
 #cpx_to_real ${procdir}/IFG/${ifg}/${ifg}.filt.diff ${procdir}/IFG/${ifg}/${ifg}.diff_mag $width 3  >> $logfile
 #cpx_to_real ${procdir}/IFG/${ifg}/${ifg}.$ifgext ${procdir}/IFG/${ifg}/${ifg}.$ifgout'_pha' $width 4  >> $logfile
 # Geocode all the data
 #geocode_back ${procdir}/IFG/${ifg}/${ifg}.$ifgext $width ${procdir}/$geodir/$master.lt_fine ${procdir}/$GEOCDIR/${ifg}/${ifg}.geo.$ifgout ${width_dem} ${length_dem} 1 1 >> $logfile
 #  geocoding complex data with Lanczos interpolation before extracting phase from it
 geocode_back ${procdir}/IFG/${ifg}/${ifg}.$ifgext $width ${procdir}/$geodir/$master.lt_fine ${procdir}/$GEOCDIR/${ifg}/${ifg}.geo.$ifgout ${width_dem} ${length_dem} 6 1 >> $logfile
 # extract only phase
 cpx_to_real ${procdir}/$GEOCDIR/${ifg}/${ifg}.geo.$ifgout ${procdir}/$GEOCDIR/${ifg}/${ifg}.geo.$ifgout'_pha' ${width_dem} 4 >> $logfile
 # convert phase to geotiff
 # data2tiff ${procdir}/$GEOCDIR/${ifg}/${ifg}.geo.$ifgout'_pha' ${width_dem} 2 ${procdir}/$GEOCDIR/${ifg}/${ifg}.geo.$ifgout'_pha.orig.tif' 0.0 >> $logfile 2>/dev/null
 data2geotiff ${procdir}/$geodir/EQA.dem_par ${procdir}/$GEOCDIR/${ifg}/${ifg}.geo.$ifgout'_pha' 2 ${procdir}/$GEOCDIR/${ifg}/${ifg}.geo.$ifgout'_pha.orig.tif' 0.0  >> $logfile 2>/dev/null
 #geocode_back ${procdir}/IFG/${ifg}/${ifg}.diff_mag $width ${procdir}/$geodir/$master.lt_fine ${procdir}/GEOC/${ifg}/${ifg}.geo.diff_mag ${width_dem} ${length_dem} 1 0 >> $logfile
 #geocode_back ${procdir}/IFG/${ifg}/${ifg}.$ifgout'_pha' $width ${procdir}/$geodir/$master.lt_fine ${procdir}/$GEOCDIR/${ifg}/${ifg}.geo.$ifgout'_pha' ${width_dem} ${length_dem} 0 0 >> $logfile
 # Convert to geotiff
 #data2geotiff ${procdir}/$geodir/EQA.dem_par ${procdir}/GEOC/${ifg}/${ifg}.geo.diff_mag 2 ${procdir}/GEOC/${ifg}/${ifg}.geo.diff_mag.tif 0.0  >> $logfile 2>/dev/null
 #data2geotiff ${procdir}/$geodir/EQA.dem_par ${procdir}/$GEOCDIR/${ifg}/${ifg}.geo.$ifgout'_pha' 2 ${procdir}/$GEOCDIR/${ifg}/${ifg}.geo.$ifgout'_pha.orig.tif' 0.0  >> $logfile 2>/dev/null
 # Some frames would need this below (quite fast):
 if [ -f $hgtfile ]; then
    gdalwarp2match.py ${procdir}/$GEOCDIR/${ifg}/${ifg}.geo.$ifgout'_pha.orig.tif' $hgtfile ${procdir}/$GEOCDIR/${ifg}/${ifg}.geo.$ifgout'_pha.orig2.tif'
    mv ${procdir}/$GEOCDIR/${ifg}/${ifg}.geo.$ifgout'_pha.orig2.tif' ${procdir}/$GEOCDIR/${ifg}/${ifg}.geo.$ifgout'_pha.orig.tif'
 fi
 # Compress
 gdal_translate -of GTiff -ot Float32 -co COMPRESS=DEFLATE -co PREDICTOR=3 ${procdir}/$GEOCDIR/${ifg}/${ifg}.geo.$ifgout'_pha.orig.tif' ${procdir}/$GEOCDIR/${ifg}/${ifg}.geo.$ifgout'_pha.tif' >> $logfile 2>/dev/null
 rm ${procdir}/$GEOCDIR/${ifg}/${ifg}.geo.$ifgout'_pha.orig.tif'
 # Create bmps
 ifg_bmp=${procdir}/$GEOCDIR/${ifg}/${ifg}.geo.$ifgout.png
 tmpifgbmp=${procdir}/$GEOCDIR/${ifg}/${ifg}.geo.$ifgout.tmp.png
 gmt grdimage ${procdir}/$GEOCDIR/${ifg}/${ifg}.geo.$ifgout'_pha.tif' -C$LiCSARpath/misc/pha.cpt -JM1 -nn+t0.1 -Q -A$ifg_bmp
 if [ $FULL -eq 1 ]; then
   convert $ifg_bmp -transparent black {procdir}/$GEOCDIR/${ifg}/${ifg}.geo.$ifgout.full.png
 fi
 #to flatten it (and fix transparency...sometimes needed..):
 convert $ifg_bmp -transparent black -resize $RESIZE'%' PNG8:$tmpifgbmp
 mv $tmpifgbmp $ifg_bmp
 
 #rasmph_pwr ${procdir}/$GEOCDIR/${ifg}/${ifg}.geo.$ifgout ${procdir}/$geodir/EQA.${master}.slc.mli ${width_dem} - - - $reducfac_dem $reducfac_dem - - - ${procdir}/$GEOCDIR/${ifg}/${ifg}.geo.$ifgout'_blk.bmp' >> $logfile
 #convert -transparent black -resize $RESIZE'%' ${procdir}/$GEOCDIR/${ifg}/${ifg}.geo.$ifgout'_blk.bmp' ${procdir}/$GEOCDIR/${ifg}/${ifg}.geo.$ifgout.png
 #if [ $FULL -eq 1 ]; then
 # gmt grdimage ${procdir}/$GEOCDIR/${ifg}/${ifg}.geo.$ifgout'_pha.tif' -C$LiCSARpath/misc/pha.cpt -JM1 -nn+t0.1 -A${procdir}/$GEOCDIR/${ifg}/${ifg}.geo.$ifgout.full.png
 # #convert -transparent black ${procdir}/$GEOCDIR/${ifg}/${ifg}.geo.$ifgout'_blk.png' ${procdir}/$GEOCDIR/${ifg}/${ifg}.geo.$ifgout.full.png
 # #rm ${procdir}/$GEOCDIR/${ifg}/${ifg}.geo.$ifgout'_blk.png'
 #fi
 #rasmph_pwr ${procdir}/GEOC/${ifg}/${ifg}.geo.diff_mag ${procdir}/$geodir/EQA.${master}.slc.mli ${width_dem} - - - $reducfac_dem $reducfac_dem - - - ${procdir}/GEOC/${ifg}/${ifg}.geo.diff_mag.bmp >> $logfile
 #rasmph_pwr ${procdir}/GEOC/${ifg}/${ifg}.geo.diff_pha ${procdir}/$geodir/EQA.${master}.slc.mli ${width_dem} - - - $reducfac_dem $reducfac_dem - - - ${procdir}/GEOC/${ifg}/${ifg}.geo.diff_pha.bmp >> $logfile

 #This was to generate only mag and pha bmp..
 #raspwr ${procdir}/GEOC/${ifg}/${ifg}.geo.diff_mag ${width_dem} - - $reducfac_dem $reducfac_dem - - - ${procdir}/GEOC/${ifg}/${ifg}.geo.diff_mag.bmp >> $logfile
 #ras_linear ${procdir}/GEOC/${ifg}/${ifg}.geo.diff_pha ${width_dem} - - $reducfac_dem $reducfac_dem - - - ${procdir}/GEOC/${ifg}/${ifg}.geo.diff_pha.bmp >> $logfile
 # Clean
 rm ${procdir}/IFG/${ifg}/${ifg}.$ifgout'_mag' ${procdir}/IFG/${ifg}/${ifg}.$ifgout'_pha' ${procdir}/$GEOCDIR/${ifg}/${ifg}.geo.$ifgout'_blk.bmp' ${procdir}/$GEOCDIR/${ifg}/${ifg}.geo.$ifgout'_pha' 2>/dev/null
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
 if [ $magcc == 1 ]; then
   cpx_to_real ${procdir}/IFG/${ifg}/${ifg}.$ifgext ${procdir}/IFG/${ifg}/${ifg}.diff_mag $width 3  >> $logfile
   geocode_back ${procdir}/IFG/${ifg}/${ifg}.diff_mag $width ${procdir}/$geodir/$master.lt_fine ${procdir}/$GEOCDIR/${ifg}/${ifg}.geo.diff_mag ${width_dem} ${length_dem} 1 0 >> $logfile
  #geocode_back ${procdir}/IFG/${ifg}/${ifg}.$ifgout'_pha' $width ${procdir}/$geodir/$master.lt_fine ${procdir}/$GEOCDIR/${ifg}/${ifg}.geo.$ifgout'_pha' ${width_dem} ${length_dem} 0 0 >> $logfile
   # Convert to geotiff
   data2geotiff ${procdir}/$geodir/EQA.dem_par ${procdir}/$GEOCDIR/${ifg}/${ifg}.geo.diff_mag 2 ${procdir}/$GEOCDIR/${ifg}/${ifg}.geo.magcc.tif 0.0  >> $logfile 2>/dev/null
   rm ${procdir}/$GEOCDIR/${ifg}/${ifg}.geo.diff_mag
   ###
   if [ -f $hgtfile ]; then
    gdalwarp2match.py ${procdir}/$GEOCDIR/${ifg}/${ifg}.geo.magcc.tif $hgtfile ${procdir}/$GEOCDIR/${ifg}/${ifg}.geo.magcc2.tif
    mv ${procdir}/$GEOCDIR/${ifg}/${ifg}.geo.magcc2.tif ${procdir}/$GEOCDIR/${ifg}/${ifg}.geo.magcc.tif
   fi
  #for compression types differences, check e.g. https://kokoalberti.com/articles/geotiff-compression-optimization-guide/
  gdal_translate -of GTiff -ot Byte -scale 0 1 0 255 -co COMPRESS=DEFLATE -co PREDICTOR=2 ${procdir}/$GEOCDIR/${ifg}/${ifg}.geo.magcc.tif ${procdir}/$GEOCDIR/${ifg}/${ifg}.geo.mag_cc.tif >> $logfile 2>/dev/null
  rm ${procdir}/$GEOCDIR/${ifg}/${ifg}.geo.magcc.tif
  ###
 fi
   if [ -f $hgtfile ]; then
    gdalwarp2match.py ${procdir}/$GEOCDIR/${ifg}/${ifg}.geo.$ifgout'_pha.orig.tif' $hgtfile ${procdir}/$GEOCDIR/${ifg}/${ifg}.geo.$ifgout'_pha.orig2.tif'
    mv ${procdir}/$GEOCDIR/${ifg}/${ifg}.geo.$ifgout'_pha.orig2.tif' ${procdir}/$GEOCDIR/${ifg}/${ifg}.geo.$ifgout'_pha.orig.tif'
   fi
 gdal_translate -of GTiff -ot Float32 -co COMPRESS=DEFLATE -co PREDICTOR=3 ${procdir}/$GEOCDIR/${ifg}/${ifg}.geo.$ifgout'_pha.orig.tif' ${procdir}/$GEOCDIR/${ifg}/${ifg}.geo.$ifgout'_pha.tif' >> $logfile 2>/dev/null
 rm ${procdir}/$GEOCDIR/${ifg}/${ifg}.geo.$ifgout'_pha.orig.tif'
 #rasmph_pwr ${procdir}/$GEOCDIR/${ifg}/${ifg}.geo.$ifgout ${procdir}/$geodir/EQA.${master}.slc.mli ${width_dem} - - - $reducfac_dem $reducfac_dem - - - ${procdir}/$GEOCDIR/${ifg}/${ifg}.geo.$ifgout'_blk.bmp' >> $logfile
 #convert -transparent black -resize $RESIZE'%' ${procdir}/$GEOCDIR/${ifg}/${ifg}.geo.$ifgout'_blk.bmp' ${procdir}/$GEOCDIR/${ifg}/${ifg}.geo.$ifgout.png
 # Create bmps
 ifg_bmp=${procdir}/$GEOCDIR/${ifg}/${ifg}.geo.$ifgout.png
 tmpifgbmp=${procdir}/$GEOCDIR/${ifg}/${ifg}.geo.$ifgout.tmp.png
 gmt grdimage ${procdir}/$GEOCDIR/${ifg}/${ifg}.geo.$ifgout'_pha.tif' -C$LiCSARpath/misc/pha.cpt -JM1 -nn+t0.1 -Q -A$ifg_bmp
 #to flatten it (and fix transparency...sometimes needed..):
 convert $ifg_bmp -transparent black -resize $RESIZE'%' PNG8:$tmpifgbmp
 mv $tmpifgbmp $ifg_bmp
 if [ $FULL -eq 1 ]; then
  gmt grdimage ${procdir}/$GEOCDIR/${ifg}/${ifg}.geo.$ifgout'_pha.tif' -C$LiCSARpath/misc/pha.cpt -JM1 -A${procdir}/$GEOCDIR/${ifg}/${ifg}.geo.$ifgout.full.png
  #convert -transparent black ${procdir}/$GEOCDIR/${ifg}/${ifg}.geo.$ifgout'_blk.png' ${procdir}/$GEOCDIR/${ifg}/${ifg}.geo.$ifgout.full.png
  #rm ${procdir}/$GEOCDIR/${ifg}/${ifg}.geo.$ifgout'_blk.png'
 fi
 rm ${procdir}/IFG/${ifg}/${ifg}.$ifgout'_mag' ${procdir}/IFG/${ifg}/${ifg}.$ifgout'_pha' ${procdir}/$GEOCDIR/${ifg}/${ifg}.geo.$ifgout'_blk.bmp' 2>/dev/null
fi


if [ $COHDO -eq 1 ]; then
if [ ! -d "${procdir}/$GEOCDIR/${ifg}" ]; then mkdir ${procdir}/$GEOCDIR/${ifg}; fi
# Unfiltered coherence: .cc, filtered: .filt.cc
for cc in cc filt.cc; do
if [ -e ${procdir}/IFG/${ifg}/${ifg}.$cc ] && [ ! -e ${procdir}/$GEOCDIR/${ifg}/${ifg}.geo.$cc.tif ]; then
  echo "Creating "$cc" coherence tiff"
  # Geocode
  geocode_back ${procdir}/IFG/${ifg}/${ifg}.$cc $width ${procdir}/$geodir/$master.lt_fine ${procdir}/$GEOCDIR/${ifg}/${ifg}.geo.$cc ${width_dem} ${length_dem} 1 0 >> $logfile
  # Convert to geotiff
  data2geotiff ${procdir}/$geodir/EQA.dem_par ${procdir}/$GEOCDIR/${ifg}/${ifg}.geo.$cc 2 ${procdir}/$GEOCDIR/${ifg}/${ifg}.geo.$cc.orig.tif 0.0  >> $logfile 2>/dev/null
   if [ -f $hgtfile ]; then
    gdalwarp2match.py ${procdir}/$GEOCDIR/${ifg}/${ifg}.geo.$cc.orig.tif $hgtfile ${procdir}/$GEOCDIR/${ifg}/${ifg}.geo.$cc.orig2.tif
    mv ${procdir}/$GEOCDIR/${ifg}/${ifg}.geo.$cc.orig2.tif ${procdir}/$GEOCDIR/${ifg}/${ifg}.geo.$cc.orig.tif
   fi
  #for compression types differences, check e.g. https://kokoalberti.com/articles/geotiff-compression-optimization-guide/
  gdal_translate -of GTiff -ot Byte -scale 0 1 0 255 -co COMPRESS=DEFLATE -co PREDICTOR=2 ${procdir}/$GEOCDIR/${ifg}/${ifg}.geo.$cc.orig.tif ${procdir}/$GEOCDIR/${ifg}/${ifg}.geo.$cc.tif >> $logfile 2>/dev/null
  rm ${procdir}/$GEOCDIR/${ifg}/${ifg}.geo.$cc.orig.tif
  # create bmps

   gmt makecpt -Cgray -T0/255/1 >${procdir}/$GEOCDIR/${ifg}/cc.cpt
   gmt grdimage ${procdir}/$GEOCDIR/${ifg}/${ifg}.geo.$cc.tif -C${procdir}/$GEOCDIR/${ifg}/cc.cpt -JM1 -nn+t0.1 -A${procdir}/$GEOCDIR/${ifg}/bbb.png
   coh_bmp=${procdir}/$GEOCDIR/${ifg}/${ifg}.geo.$cc.png
   convert -transparent black -resize $RESIZE'%' ${procdir}/$GEOCDIR/${ifg}/bbb.png PNG8:$coh_bmp
   
  #new version of gamma shows coherence in colour... using old-school cpxfiddle as workaround
  #rascc ${procdir}/GEOC/${ifg}/${ifg}.geo.cc - ${width_dem} - - - $reducfac_dem $reducfac_dem 0 1 - - - ${procdir}/GEOC/${ifg}/${ifg}.geo.cc_blk.bmp >> $logfile
 # module load doris
  #byteswap -o ${ifg}.geo.cc.tmp 4 ${procdir}/$GEOCDIR/${ifg}/${ifg}.geo.cc >/dev/null
  #will use the gamma swapper here..
  #swap_bytes ${procdir}/$GEOCDIR/${ifg}/${ifg}.geo.cc ${ifg}.geo.cc.tmp 4 >/dev/null
  #cpxfiddle -q normal -w ${width_dem} -f r4 -o sunraster -c gray -M $reducfac_dem/$reducfac_dem -r 0/0.9 ${ifg}.geo.cc.tmp > ${procdir}/$GEOCDIR/${ifg}/${ifg}.geo.cc_blk.ras 2>/dev/null
  #convert ${procdir}/$GEOCDIR/${ifg}/${ifg}.geo.cc_blk.ras ${procdir}/$GEOCDIR/${ifg}/${ifg}.geo.cc_blk.bmp
  rm -f ${procdir}/$GEOCDIR/${ifg}/${ifg}.geo.cc_blk.ras ${ifg}.geo.$cc.tmp 2>/dev/null
  # Need to remove the black border, but the command below is no good as it removes the black parts of the coherence!
  #convert -transparent black -resize 30% ${procdir}/GEOC/${ifg}/${ifg}.geo.cc_blk.bmp ${procdir}/GEOC/${ifg}/${ifg}.geo.cc.bmp
  #convert -transparent black -resize $RESIZE'%' ${procdir}/$GEOCDIR/${ifg}/${ifg}.geo.cc_blk.bmp ${procdir}/$GEOCDIR/${ifg}/${ifg}.geo.cc.png
  if [ $FULL -eq 1 ]; then
   convert -transparent black ${procdir}/$GEOCDIR/${ifg}/bbb.png ${procdir}/$GEOCDIR/${ifg}/${ifg}.geo.$cc.full.png
  fi
  rm ${procdir}/$GEOCDIR/${ifg}/${ifg}.geo.cc_blk.bmp 2>/dev/null
  rm ${procdir}/$GEOCDIR/${ifg}/bbb.png ${procdir}/$GEOCDIR/${ifg}/cc.cpt ${procdir}/$GEOCDIR/${ifg}/${ifg}.geo.$cc 2>/dev/null
fi
done
fi

if [ $MLIDO -eq 1 ]; then
#finally amplitudes/intensities
mkdir $procdir/$GEOCDIR.MLI 2>/dev/null
for im in `echo $ifg | sed 's/_/ /'`; do
if [ -e ${procdir}/RSLC/$im/$im.rslc.mli ] && [ ! -d ${procdir}/$GEOCDIR.MLI/$im ]; then
 mkdir -p ${procdir}/$GEOCDIR.MLI/$im
 echo "preparing GeoTIFF and PNG for intensity image of "$im
 #geocode MLI
 #echo geocode_back ${procdir}/RSLC/$im/$im.rslc.mli $width ${procdir}/$geodir/$master.lt_fine ${procdir}/$GEOCDIR.MLI/$im/$im.geo.mli ${width_dem} ${length_dem} 0 0 #>> $logfile
 geocode_back ${procdir}/RSLC/$im/$im.rslc.mli $width ${procdir}/$geodir/$master.lt_fine ${procdir}/$GEOCDIR.MLI/$im/$im.geo.mli ${width_dem} ${length_dem} 0 0 >> $logfile
 #convert MLI to geotiff
 #echo data2geotiff ${procdir}/$geodir/EQA.dem_par ${procdir}/$GEOCDIR.MLI/$im/$im.geo.mli 2 ${procdir}/$GEOCDIR.MLI/$im/$im.geo.mli.tif 0.0  #>> $logfile 2>/dev/null
 data2geotiff ${procdir}/$geodir/EQA.dem_par ${procdir}/$GEOCDIR.MLI/$im/$im.geo.mli 2 ${procdir}/$GEOCDIR.MLI/$im/$im.geo.mli.orig.tif 0.0  >> $logfile 2>/dev/null
   if [ -f $hgtfile ]; then
    gdalwarp2match.py ${procdir}/$GEOCDIR.MLI/$im/$im.geo.mli.orig.tif $hgtfile ${procdir}/$GEOCDIR.MLI/$im/$im.geo.mli.orig2.tif
    mv ${procdir}/$GEOCDIR.MLI/$im/$im.geo.mli.orig2.tif ${procdir}/$GEOCDIR.MLI/$im/$im.geo.mli.orig.tif
   fi
  #for compression types differences, check e.g. https://kokoalberti.com/articles/geotiff-compression-optimization-guide/
  gdal_translate -of GTiff -co COMPRESS=DEFLATE -co PREDICTOR=3 ${procdir}/$GEOCDIR.MLI/$im/$im.geo.mli.orig.tif ${procdir}/$GEOCDIR.MLI/$im/$im.geo.mli.tif >> $logfile 2>/dev/null
  
 #generate raster preview
 raspwr ${procdir}/$GEOCDIR.MLI/$im/$im.geo.mli ${width_dem} - - $reducfac_dem $reducfac_dem - - - ${procdir}/$GEOCDIR.MLI/$im/$im.geo.mli.bmp 0 - >> $logfile
 convert -transparent black -resize $RESIZE'%' ${procdir}/$GEOCDIR.MLI/$im/$im.geo.mli.bmp ${procdir}/$GEOCDIR.MLI/$im/$im.geo.mli.png
 rm ${procdir}/$GEOCDIR.MLI/$im/$im.geo.mli.bmp ${procdir}/$GEOCDIR.MLI/$im/$im.geo.mli  ${procdir}/$GEOCDIR.MLI/$im/$im.geo.mli.orig.tif 2>/dev/null
fi 
done
fi

rm gmt.history 2>/dev/null

echo "done"


exit


if [ $mask -eq 1 ]; then
  #exit
  frame=`basename $procdir`
  tr=`echo $frame | sed 's/^0//' | sed 's/^0//'`
  if [ -f $LiCSAR_public/$tr/$frame/metadata/$frame.geo.landmask.tif ]; then
   maskfile=$LiCSAR_public/$tr/$frame/metadata/$frame.geo.landmask.tif
  elif [ -f $procdir/$frame.geo.landmask.tif ]; then
   maskfile=$procdir/$frame.geo.landmask.tif
  else 
   wget --no-check-certificate -O $procdir/$frame.geo.landmask.tif https://gws-access.ceda.ac.uk/public/nceo_geohazards/LiCSAR_products/$tr/$frame/metadata/$frame.geo.landmask.tif 2>/dev/null
   if [ -f $procdir/$frame.geo.landmask.tif ]; then maskfile=$procdir/$frame.geo.landmask.tif; else maskfile=''; fi
  fi
  if [ ! -z $maskfile ]; then
   echo "applying mask"
   for x in `ls $procdir/GEOC/${ifg}/*tif`; do
    mv $x $x.temp.tif
    gdal_calc.py -A $x.temp.tif -B $maskfile --NoDataValue=0 --outfile=$x --calc="A*B" >/dev/null 2>/dev/null
    rm $x.temp.tif
   done
  fi
fi
#Move it to the public area(?)
#for filename in .......; do
# if [ ! -f ${publicdir}/ ]

