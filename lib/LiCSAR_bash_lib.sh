#!/bin/bash
#echo "LiCSAR Bash Libraries 2016 are loaded"
#########################################################
# LICSAR Bash Shell Function Library
#########################################################
#
# Author: Pablo J. Gonzalez <p.j.gonzalez@leeds.ac.uk>
# Contributions by: Emma Hatton and Richard Walters, and Milan Lazecky
#
# Copyright 2015-2020
# Released under the current GPL version.
#
# Description:
#   This is a shell script library. It contains functions 
#   that can be called by programs that include (source)
#    the LiCSAR InSAR software system. 
# Updates:
#    01/2021: add functions to do solid Earth tides correction
# By simply sourcing this library, you can use all available functions as 
# documented on the project page [to be done!].
#
##########################################################

#fix for bsub -> slurm
a=`which bsub 2>/dev/null`
if [ -z $a ]; then
 echo "setting bsub2slurm"
 alias bsub=bsub2slurm.sh
fi


function slurm_get_running_jobs() {
  for jobid in `bjobs | grep ' R ' | gawk {'print $1'}`; do
   slurm_get_jobname $jobid
  done
}

function slurm_get_jobname() {
    jobid=$1
    scontrol show jobid $jobid | grep JobName | cut -d '=' -f3
}

function touchscratch() {
    dir=$1
    if [ `echo $dir | grep -c scratch` == 0 ]; then
      echo "Usage: touchscratch \$BATCH_CACHE_DIR";
      echo "(or another directory on scratch disk)"
      return 0;
    fi;
    echo "touching contents of directory "$dir
    for x in `find $dir`; do touch $x; done
    echo "done"
}

function clean_public() {
    cd $LiCSAR_public
    for track in `seq 1 175`; do echo $track; 
      for frame in `ls $track | grep '_'`; do
         if [ `ls $track/$frame | wc -l` -lt 2 ]; then
            echo "frame "$frame" seems wrongly intialized, see list:"
            ls $track/$frame
         fi
         if [ -d $track/$frame/interferograms/geo ]; then rm -rf $track/$frame/interferograms/geo; fi
         for ifg in `ls $track/$frame/interferograms`; do
            if [ `ls $track/$frame/interferograms/$ifg | wc -l` -lt 2 ]; then
               echo "Removing empty folder of"
               echo $track/$frame/interferograms/$ifg
               rm -rf $track/$frame/interferograms/$ifg
            fi
         done
      done
    done
}

function create_LOS_tide() {
  if [ "$#" == "2" ]; then
    local frame=$1
    local date1=$2
    #fix for input with yyyymmdd
    datem1=`date -d $date1 +%Y-%m-%d`
    track=`echo $frame | cut -c -3 | sed 's/^0//' | sed 's/^0//'`
    #get center time
    epochpath=$LiCSAR_public/$track/$frame/epochs/$date1
    if [ -f $epochpath/$date1'.tide.geo.tif' ]; then
      echo "tide file for "$date1" already exists"
      return 0
    fi
    echo "generating solid Earth tide in LOS for date "$date1
    source $LiCSAR_public/$track/$frame/metadata/metadata.txt
    U=$LiCSAR_public/$track/$frame/metadata/$frame.geo.U.tif
    E=$LiCSAR_public/$track/$frame/metadata/$frame.geo.E.tif
    N=$LiCSAR_public/$track/$frame/metadata/$frame.geo.N.tif
    mkdir -p $epochpath
    gmt earthtide -T$datem1'T'$center_time -G$epochpath/tmp_tides.%s.nc -R$U -Cv,e,n 2>/dev/null
    gmt grdmath $U $epochpath/tmp_tides.v.nc MUL $E $epochpath/tmp_tides.e.nc MUL ADD $N $epochpath/tmp_tides.n.nc MUL ADD = $epochpath/$date1'.tide.geo.tif'=gd:GTiff
    rm $epochpath/tmp_tides.*.nc
  else
    echo "Usage: create_LOS_tide frame_id date"
    echo "(date must be yyyymmdd)"
    return 0
  fi
};


function prepare_landmask() {
    infile=$1
    frame=$2
    filedir=`dirname $infile`
    tr=`track_from_frame $frame`
    if [ -f $LiCSAR_public/$tr/$frame/metadata/$frame.geo.landmask.tif ]; then
        maskfile=$LiCSAR_public/$tr/$frame/metadata/$frame.geo.landmask.tif
    elif [ -f $frame.geo.landmask.tif ]; then
        maskfile=$frame.geo.landmask.tif
    elif [ -f $filedir/$frame.geo.landmask.tif ]; then
        maskfile=$filedir/$frame.geo.landmask.tif
    else 
        wget --no-check-certificate -O $filedir/$frame.geo.landmask.tif https://gws-access.jasmin.ac.uk/public/nceo_geohazards/LiCSAR_products/$tr/$frame/metadata/$frame.geo.landmask.tif 2>/dev/null
        if [ -f $filedir/$frame.geo.landmask.tif ]; then
            if [ `ls -l $filedir/$frame.geo.landmask.tif | gawk {'print $5'}` -eq 0 ]; then
                rm $filedir/$frame.geo.landmask.tif; maskfile='';
            else
                maskfile=$filedir/$frame.geo.landmask.tif;
            fi;
        else maskfile=''; fi
    fi
    landmask=0
    if [ ! -z $maskfile ]; then
        maskmin=`gdalinfo -stats $maskfile 2>/dev/null | grep ICS_MINIMUM | cut -d '=' -f2`
        if [ $maskmin -eq 0 ]; then
            #resample to the AOI (fix for the missing bursts etc.)
            tmpfile1=$infile.temp.nc
            gmt grdedit -T $infile -G$tmpfile1
            gmt grdedit -T -R$tmpfile1 $maskfile -G$filedir/tempmask.nc
            #gmt grdsample $maskfile -Gtempmask.nc -R$unw_tif -nn+t0.1 2>/dev/null
            #not working well! skipping
            #gmt grdconvert $maskfile $filedir/tempmask.nc
            #gmt grdcut -N -R$unw_tif -Gtempmask.nc $maskfile 2>/dev/null
            if [ -f $filedir/tempmask.nc ]; then
                #echo "will apply landmask"
                masknc=$filedir/tempmask.nc 
                landmask=1
            fi
            
        fi
    fi
    if [ $landmask -eq 1 ]; then
      gmt grdmath $tmpfile1 $masknc MUL 0 NAN = $filedir/tempmasked.nc
      ls $filedir/tempmasked.nc
      rm $masknc 2>/dev/null
    fi
    if [ $landmask -eq 2 ]; then
        #should apply landmask  -  file $maskfile
        rm $filedir/tempmasked.nc 2>/dev/null
        gmt grdclip $infile -G$filedir/tempinfile1.nc -SrNaN/0
        #gmt grdsample tempinfile1.nc -Gtempinfile2.nc -R$masknc -nn+t0.1 2>/dev/null
        #gmt grdcut -N1 tempinfile1.nc -Gtempinfile2.nc -R$infile
        #just convert mask to pixel registration..
        gmt grdedit $masknc -G$masknc.nc -T -R$filedir/tempinfile1.nc
        mv $masknc.nc $masknc
        gmt grdcut -N1 $masknc -G$filedir/tempmask2.nc -R$filedir/tempinfile1.nc #$infile
        if [ ! `gmt grdinfo $filedir/tempmask2.nc | grep "node registration used" | gawk {'print $2'}` == 'Pixel' ]; then
         #just convert to pixel registration..
         gmt grdedit $filedir/tempinfile1.nc -G$filedir/tempinfile1.nc.nc -T -R$filedir/tempinfile1.nc
         mv $filedir/tempinfile1.nc.nc $filedir/tempmask2.nc
        fi
        gmt grdmath -N $filedir/tempinfile1.nc $filedir/tempmask2.nc MUL = $filedir/temp2.nc
        #gmt grdmath -N tempinfile2.nc $masknc MUL = temp2.nc
        gmt grdclip $filedir/temp2.nc -G$filedir/tempmasked.nc -Sr0/NaN
        ls $filedir/tempmasked.nc
        rm $filedir/temp2.nc $filedir/tempinfile1.nc $filedir/tempmask2.nc $filedir/tempmask.nc 2>/dev/null
    fi
    if [ ! -z $tmpfile1 ]; then
       rm $tmpfile1 2>/dev/null
    fi
}

function prepare_geo_56m() {
   if [ ! -z $1 ]; then
    framepath=$1
    frame=`basename $framepath`
    #this will work for outres=0.0005
    #e.g. framepath=/home/home02/earmla/licsar/eq/frames/$frame
    demdir=$framepath/DEM
    geodir=$framepath/geo
    master=`ls $geodir/20??????.hgt | rev | cut -d '/' -f1 | rev | cut -d '.' -f1`
    masterdir=$framepath/RSLC/$master
    mv $geodir $geodir'.backup_lowres'
    mkdir $geodir
    python3 -c "import datetime as dt; from LiCSAR_lib.coreg_lib import geocode_dem; a='"$master"'; masterdate = dt.date(int(a[:4]),int(a[4:6]),int(a[6:8])); \
geocode_dem('"$masterdir"','"$geodir"','"$demdir"','"$procdir"',masterdate,0.0005)"
   else
    echo "Usage: prepare_geo_56m full_path_to_framedir"
    echo "e.g. $LiCSAR_public/51/051x_..."
    echo "this function will just generate geo folder that will have 56m resolution lookup tables for geocoding to that resolution.."
    return 0
  fi
}

 
function prepare_hillshade() {
    infile=$1
    filedir=`dirname $infile`
    frame=$2
    tr=`track_from_frame $frame`
    hgtfile=$LiCSAR_public/$tr/$frame/metadata/$frame.geo.hgt.tif
hillshadefile=$LiCSAR_public/$tr/$frame/metadata/$frame.geo.hillshade.nc
if [ ! -f $hgtfile ]; then
  mkdir -p ../../geo
  hgtfile=../../geo/$frame.geo.hgt.tif
  if [ ! -f $hgtfile ]; then
   #try for the last time, download it..
   wget --no-check-certificate -O $hgtfile https://gws-access.ceda.ac.uk/public/nceo_geohazards/LiCSAR_products/$tr/$frame/metadata/$frame.geo.hgt.tif 2>/dev/null
  fi
  #if [ ! -f $hgtfile ]; then
  # echo "warning, no hgt tiff file found. there will be no hillshade"
  #fi
fi

if [ ! -f $hillshadefile ]; then
 if [ ! -d $LiCSAR_public/$tr/$frame/metadata ]; then
  #in such case do it only locally
  hillshadefile=$filedir/hillshade.nc
 fi
fi

if [ ! -f $hillshadefile ]; then
   #do hillshade to unw..
   #echo "(generating hillshade)"
   if [ -f $hgtfile ]; then
    opass=`echo $frame | cut -c 4`
    if [ $opass == 'A' ]; then
      deg=258
    else
      deg=102
    fi
    #gmt grdsample $hgtfile -Gtemphgt.nc -R$unw_tif 2>/dev/null
    gmt grdgradient -A$deg -Nt1 -G$hillshadefile $hgtfile
   fi
fi

if [ -f $hillshadefile ]; then
      hillshadenc=$filedir/temphillshade.nc
      rm $hillshadenc 2>/dev/null
      gmt grdsample $hillshadefile -G$hillshadenc.tmp -R$infile -nn+t0.1 2>/dev/null
      gmt grdcut -N1 $hillshadenc.tmp -G$hillshadenc -R$infile
      rm $hillshadenc.tmp
  fi
  if [ -f $hillshadenc ]; then
    echo $hillshadenc
  fi
}


function create_colourbar_unw() {
    infile=$1
  minmaxcolour=`gmt grdinfo -T+a1+s $infile`
  #create legend
  #need to prepare a colorbar based on these values!!!!
  #you know... the rounding here is not really important... just a colourbar.. or not?
  mincol=`echo $minmaxcolour | cut -d '/' -f1 | cut -d 'T' -f2`
  maxcol=`echo $minmaxcolour | cut -d '/' -f2 `
  #expecting sentinel
  minval=`python -c 'print(round('$mincol'*5.546/(4*3.14159265)))'`
  maxval=`python -c 'print(round('$maxcol'*5.546/(4*3.14159265)))'`
  #add also real min and max values
  minmaxreal=`gmt grdinfo -T $infile`
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
  convert -font helvetica -fill black -pointsize 40 -draw "text "$xsize",115 '"$minval"'" $LiCSARpath/misc/scalebar_unwrapped_empty.png temp_scale_unw.png
  convert -font helvetica -fill black -pointsize 40 -draw "text 1100,115 '"$maxval" cm'" temp_scale_unw.png scalebar_unwrapped.png
  mv scalebar_unwrapped.png temp_scale_unw.png
  #add real values
  convert -font helvetica -fill black -pointsize 35 -draw "text "$xsize",165 '[min "$minrealval" cm]'" temp_scale_unw.png scalebar_unwrapped.png
  mv scalebar_unwrapped.png temp_scale_unw.png
  convert -font helvetica -fill black -pointsize 35 -draw "text 1020,165 '[max "$maxrealval" cm]'" temp_scale_unw.png scalebar_unwrapped.png
  rm temp_scale_unw.png
  
  scalebarfile='scalebar_unwrapped.png'
  echo $scalebarfile
}

function create_colourbar_m() {
    infile=$1
    #code=rngoff or azioff
    code=$2
    cutoffperc=$3 # a=cutoff percentage for values - 25% can be ok
  minmaxcolour=`gmt grdinfo -T+a$cutoffperc'+s' $infile`  
  #create legend
  #need to prepare a colorbar based on these values!!!!
  #you know... the rounding here is not really important... just a colourbar.. or not?
  minval=`echo $minmaxcolour | cut -d '/' -f1 | cut -d 'T' -f2`
  maxval=`echo $minmaxcolour | cut -d '/' -f2 `
  
  minval=`python -c 'print(round('$minval'))'` #'*5.546/(4*3.14159265)))'`
  maxval=`python -c 'print(round('$maxval'))'` #'*5.546/(4*3.14159265)))'`
  #add also real min and max values
  minmaxreal=`gmt grdinfo -T $infile`
  minrealval=`echo $minmaxreal | cut -d 'T' -f2 | cut -d '/' -f1 | cut -d '.' -f1`
  maxrealval=`echo $minmaxreal | cut -d '/' -f2 | cut -d '.' -f1`
  #minrealval=`python -c 'print(round('$minreal'*5.546/(4*3.14159265)))'`
  #maxrealval=`python -c 'print(round('$maxreal'*5.546/(4*3.14159265)))'`
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
  convert -font helvetica -fill black -pointsize 40 -draw "text "$xsize",115 '"$minval"'" $LiCSARpath/misc/scalebar_$code'off_empty'.png $infile.temp_scale.png
  convert -font helvetica -fill black -pointsize 40 -draw "text 1100,115 '"$maxval" m'" $infile.temp_scale.png $infile.scalebar_m.png
  mv $infile.scalebar_m.png $infile.temp_scale.png
  #add real values
  convert -font helvetica -fill black -pointsize 35 -draw "text "$xsize",165 '[min "$minrealval" m]'" $infile.temp_scale.png $infile.scalebar_m.png
  mv $infile.scalebar_m.png $infile.temp_scale.png
  convert -font helvetica -fill black -pointsize 35 -draw "text 1020,165 '[max "$maxrealval" m]'" $infile.temp_scale.png $infile.scalebar_m.png
  rm $infile.temp_scale.png
  
  scalebarfile=`dirname $infile`/'scalebar_'$code'.png'
  mv $infile.scalebar_m.png $scalebarfile
  echo $scalebarfile
}



function create_preview_offsets() {
    # usage: infile (either geo.rng.tif or geo.azi.tif) and frame
    # if third parameter is a number, it will be used as a cutoff + will not delete scale bar (for kmz use)
  if [ ! -z $1 ]; then
    local infile=$1
    code=`basename $infile | cut -d '.' -f3`
    echo "generating preview for "$infile
    outfile=`echo $infile | rev | cut -c 4- | rev`png
    #extracmd_convert='-resize 30%'
    extracmd=''
    extracmd_convert=''
    #if [ ! -z $2 ]; then
    # if [ $2 == 0 ]; then
    #  extracmd_convert='-resize 30%'
    # fi
    #fi
   if [ ! -f `basename $outfile .png`.masked.tif ]; then
    if [ ! -z $2 ]; then
       frame=$2
       tr=`track_from_frame $frame`
       #echo "trying mask and include hillshade"
       maskedfile=`prepare_landmask $infile $frame`
       if [ ! -z $maskedfile ]; then
        infile=$maskedfile
       fi
    fi
   else
    infile=`basename $outfile .png`.masked.tif
   fi
    if [ ! -z $3 ]; then cutoff=$3; else cutoff=20; fi
    minoff=-10
    maxoff=10
    echo "limiting offsets to "$minoff"/"$maxoff" m, in preview cutting "$cutoff" percent of tails"
    gmt grdclip $infile -G`basename $outfile .png`.masked.tif=gd:Gtiff -Sr0/NaN -Sb$minoff/NaN -Sa$maxoff/NaN
    #gmt grdconvert $infile.masked.nc -G$outfile.masked.tif:GTiff
    infile=`basename $outfile .png`.masked.tif
    barpng=`create_colourbar_m $infile $code $cutoff`
    minmaxcolour=`gmt grdinfo -T+a$cutoff'+s' $infile` # must remain same as in create_colourbar_m !
    gmt makecpt -C$LiCSARpath/misc/colourmap.cpt -Iz $minmaxcolour/0.025 >`dirname $outfile`/$code.cpt
    gmt grdimage $infile -C`dirname $outfile`/$code.cpt $extracmd -JM1 -Q -nn+t0.1 -A$outfile.tt.png
    #convert $extracmd_convert $outfile.tt.png PNG8:$outfile; rm $outfile.tt.png
    convert $outfile.tt.png PNG8:$outfile; rm $outfile.tt.png
   if [ ! -z $frame ]; then
    #convert $outfile.tt.png PNG8:$outfile; rm $outfile.tt.png
    #if [ ! -z $frame ]; then
    if [ `echo $frame | cut -c 4` == 'A' ]; then grav='southeast'; else grav='southwest'; fi
   else
    # no frame? firmly bottom left corner
    grav='southwest'
   fi
   convert $outfile -resize 680x \( $barpng -resize 400x  -background none -gravity center \) -gravity $grav -geometry +7+7 -composite -flatten -transparent black $outfile.temp.png
   convert $outfile.temp.png -transparent black $extracmd_convert PNG8:$outfile
   if [ ! -z $4 ]; then
     if [ $code == 'rng' ]; then
      TN="range_offsets"
     else
      TN="azimuth_offsets"
     fi
     echo 'debug: generating kml for '$code
     # overcoming the annoying gridline/pixel registration issue... finally
     gmt grdconvert $infile -G$infile.nc
     gmt grdedit -T $infile.nc
     gmt grd2kml -Ag -C`dirname $outfile`/$code.cpt -nn+t0.1 -T$TN -N$TN $extracmd $infile.nc 2>/dev/null
     rm $infile.nc
   else
    rm `dirname $outfile`/$code.cpt
    rm $barpng $outfile.temp.png
   fi
  else
    echo "Usage: create_preview_offsets ..geo.rng/azi.tif [frame] [cutoff] [tokml?]"
    echo "(can be either geotiff or nc/grd)"
    echo "cutoff is 20 percent trim by default. tokml=1 means create kml"
    return 0
  fi
}

function create_preview_vel() {
	# so the colour scale was generated as:
	# 
  if [ ! -z $1 ]; then
    local velfile=$1
    echo "generating preview for "$velfile
    outfile=`echo $velfile | rev | cut -c 4- | rev`png
    extracmd_convert=''
    if [ ! -z $2 ]; then
     if [ $2 == 0 ]; then
      extracmd_convert='-resize 30%'
     fi
    fi
    gmt grdimage $velfile -C$LiCSARpath/misc/vel.cpt -JM1 -nn+t0.1 -Q -A$outfile.temp.png
    #to flatten it (and fix transparency...sometimes needed..):
    convert $outfile.temp.png -transparent black $extracmd_convert PNG8:$outfile
    rm $outfile.temp.png
    rm gmt.history 2>/dev/null
  else
    echo "Usage: create_preview_vel velocity_tif_file [hires?]"
    echo "(can be either geotiff or nc/grd; if hires=0, it will NOT keep hires..)"
    return 0
  fi
}


function create_preview_unwrapped() {
  if [ ! -z $1 ]; then
    local unwfile=$1
    echo "generating preview for "$unwfile
    origfile=$unwfile
    #local RESIZE=30
    outfile=`echo $unwfile | rev | cut -c 4- | rev`png
    # correct for median
    unwfile=$origfile.tmp.tif
    gmt grdmath $origfile $origfile MEDIAN SUB = $unwfile=gd:Gtiff
  extracmd=''
  if [ ! -z $2 ]; then
   frame=$2
   tr=`track_from_frame $frame`
   #echo "trying mask and include hillshade"
   maskedfile=`prepare_landmask $unwfile $frame`
   if [ ! -z $maskedfile ]; then
    unwfile=$maskedfile
   fi
   hillshade=`prepare_hillshade $unwfile $frame`
   if [ ! -z $hillshade ]; then
    extracmd='-I'$hillshade
   fi
  fi
  barpng=`create_colourbar_unw $unwfile`
  minmaxcolour=`gmt grdinfo -T+a1+s $unwfile`
  gmt makecpt -C$LiCSARpath/misc/colourmap.cpt -Iz $minmaxcolour/0.025 >$outfile.unw.cpt
  if [ -z $3 ]; then
   gmt grdimage $unwfile -C$outfile.unw.cpt $extracmd -JM1 -Q -nn+t0.1 -A$outfile.tt.png
   convert $outfile.tt.png PNG8:$outfile; rm $outfile.tt.png
   if [ ! -z $frame ]; then
    if [ `echo $frame | cut -c 4` == 'A' ]; then grav='southeast'; else grav='southwest'; fi
   else
    # no frame? firmly bottom left corner
    grav='southwest'
   fi
   convert $outfile -resize 680x \( $barpng -resize 400x  -background none -gravity center \) -gravity $grav -geometry +7+7 -composite -flatten -transparent black $outfile.sm.png
   #save only the small preview..
   mv $outfile.sm.png $outfile
   rm $barpng $outfile.unw.cpt
  else
   #echo "preparing for kml"
   gmt grd2kml -Ag -C$outfile.unw.cpt -nn+t0.1 -Tunwrapped_ifg -Nunwrapped_ifg $extracmd $unwfile 2>/dev/null
  fi
  rm $maskedfile $hillshade $unwfile $outfile.unw.cpt 2>/dev/null
  rm gmt.history 2>/dev/null
  else
    echo "Usage: create_preview_unwrapped unwrapped_ifg [frame] [to kmz?]"
    echo "(can be either geotiff or nc/grd; if frame is provided, it will use mask/hillshade)"
    return 0
  fi
};


function create_landmask() {
  if [ ! -z $1 ]; then
    local tif=$1
    masktype="0/1/0/1/0"
    echo $masktype " for .. <ocean>/<land>/<lake>/<island>/<pond>"
    gmt grdlandmask -Glandmask.tif=gd:GTiff -R$tif -Df -N$masktype
  else
    echo "Usage: create_landmask tiffile"
    echo "(can be either geotiff or nc/grd)"
    return 0
  fi
}

function create_preview_coh() {
  if [ ! -z $1 ]; then
    local cohfile=$1
    extracmd_convert=''
    if [ ! -z $2 ]; then
     extracmd_convert='-resize 30%'
    fi
    outfile=`echo $cohfile | sed 's/_pha//' | rev | cut -c 4- | rev`png
    gmt grdimage $cohfile -C$LiCSARpath/misc/cc.cpt -JM1 -nn+t0.1 -Q -A$outfile.temp.png
    #to flatten it (and fix transparency...sometimes needed..):
    convert $outfile.temp.png -transparent black $extracmd_convert PNG8:$outfile
    rm $outfile.temp.png
  else
    echo "Usage: create_preview_coh cohfile [hires?]"
    echo "(can be either geotiff or nc/grd; if hires=1, it will keep hires)"
    return 0
  fi
}

function create_preview_wrapped() {
  if [ ! -z $1 ]; then
    local ifgfile=$1
    extracmd_convert=''
    if [ ! -z $2 ]; then
     if [ $2 == 0 ]; then
      extracmd_convert='-resize 30%'
     fi
    fi
    outfile=`echo $ifgfile | sed 's/_pha//' | rev | cut -c 4- | rev`png
    gmt grdimage $ifgfile -C$LiCSARpath/misc/pha.cpt -JM1 -nn+t0.1 -Q -A$outfile.temp.png
    #to flatten it (and fix transparency...sometimes needed..):
    convert $outfile.temp.png -transparent black $extracmd_convert PNG8:$outfile
    rm $outfile.temp.png
  else
    echo "Usage: create_preview_wrapped ifgfile [hires?]"
    echo "(can be either geotiff or nc/grd; if hires=1, it will keep hires)"
    return 0
  fi
};

function track_from_frame() {
 #input is $frame
 echo $1 | cut -c -3 | sed 's/^0//' | sed 's/^0//'
};


function correct_ifg_gacos_public() {
  if [ "$#" == "3" ]; then
    local frame=$1
    local ifg=$2
    local ext=$3
    track=`echo $frame | cut -c -3 | sed 's/^0//' | sed 's/^0//'`
    ifgpath=$LiCSAR_public/$track/$frame/interferograms/$ifg
    infile=$ifgpath/$ifg.geo.$ext.tif
    date1=`echo $ifg | cut -d '_' -f1`
    date2=`echo $ifg | cut -d '_' -f2`
    gacos1=$LiCSAR_public/$track/$frame/epochs/$date1/$date1'.sltd.geo.tif'
    gacos2=$LiCSAR_public/$track/$frame/epochs/$date2/$date2'.sltd.geo.tif'
    if [ ! -f $infile ]; then
      echo "interferogram does not exist"
      return 0
    fi
    for gacosf in $gacos1 $gacos2; do
     if [ ! -f $gacosf ]; then
      echo "gacos file "$gacosf" does not exist"
      return 0
     fi
    done
    #U=$LiCSAR_public/$track/$frame/metadata/$frame.geo.U.tif
    outfile=$ifgpath/$ifg.geo.$ext.gacos.tif
    #wavelength of Sentinel-1
    #wavelength=0.055465763
    #mcoef=4*PI/$wavelength
    #mcoef=226.56
    echo "correcting the interferogram to "$outfile
    if [ $ext == 'unw' ]; then
     gmt grdmath $infile'=gd:Gtiff+n0' 0 NAN $gacos2 $gacos1 SUB SUB = $outfile'=gd:Gtiff'
     create_preview_unwrapped $outfile
    else
     #sltd are already in radians
     gmt grdmath $infile'=gd:Gtiff+n0' 0 NAN $gacos2 $gacos1 SUB SUB WRAP = $outfile'=gd:Gtiff'
     #gmt grdmath $infile'=gd:Gtiff+n0' 0 NAN $gacos2 $gacos1 SUB 226.56 MUL 0 NAN $U DIV SUB WRAP = $outfile'=gd:Gtiff'
     create_preview_wrapped $outfile
    fi
  else
    echo "Usage: correct_ifg_gacos_public frame_id ifg ext"
    echo "(ifg must be yyyymmdd_yyyymmdd)"
    echo "(ext should be either diff_pha, diff_unfiltered_pha or unw)"
    return 0
  fi
};


function cdproc() {
    if [ -z $1 ]; then
        echo "please provide frame id";
        return 0
    fi
    frame=$1
    tr=`track_from_frame $frame`
    cd $LiCSAR_procdir/$tr/$frame
}

function cdpub() {
    if [ -z $1 ]; then
        echo "please provide frame id";
        return 0
    fi
    frame=$1
    tr=`track_from_frame $frame`
    cd $LiCSAR_public/$tr/$frame
}

function recreate_frame_to_hires() {
    #to be run for all tienshan frames that were initialised 'normally'
    if [ -z $1 ]; then
        echo "please provide frame id";
        return 0
    fi
    echo "WARNING - expecting this frame to be from Tien Shan, check local_config.py if this is not wanted"
    local frame=$1
    echo "regenerating the frame files to make medium resolution outputs (56 m)"
    tr=`track_from_frame $frame`
    master=`get_master $frame`
    m=$master
    fdir=$LiCSAR_procdir/$tr/$frame
    cd $fdir
    echo "outres=0.0005" > local_config.py
    echo "tienshan=1" >> local_config.py
    mv geo backup.geo.110m
    mkdir geo
    echo "master is "$master
    python3 -c "from LiCSAR_lib.coreg_lib import geocode_dem; geocode_dem('RSLC/"$master"','geo','DEM','.','"$master"',0.0005)"
    
    LiCSAR_05_mk_angles_master
    echo "Generating E-N-U files"
    submit_lookangles.py -f $frame -t $tr
    echo "Generating land mask (would remove heights below 0 m)"
    landmask=$LiCSAR_procdir/$tr/$frame/geo/landmask
    hgtgeo=$LiCSAR_public/$tr/$frame/metadata/$frame.geo.hgt.tif
    gmt grdlandmask -G$landmask.tif=gd:GTiff -R$hgtgeo -Df -N0/1/0/1/0
    cp $landmask.tif $LiCSAR_public/$tr/$frame/metadata/$frame.geo.landmask.tif
    echo "Generating master MLI geotiff"
    create_geoctiffs_to_pub.sh -M `pwd` $m
    mkdir -p $LiCSAR_public/$tr/$frame/epochs 2>/dev/null
    rm -r $LiCSAR_public/$tr/$frame/epochs/$m 2>/dev/null
    mv GEOC.MLI/$m $LiCSAR_public/$tr/$frame/epochs/.
    rmdir GEOC.MLI 2>/dev/null
    echo "Generating public metadata file"
    submit_to_public_metadata.sh $frame
    ifgdir=$LiCSAR_public/$tr/$frame/interferograms
    if [ `ls $ifgdir 2>/dev/null | wc -l` -gt 0 ]; then
        echo "backing up ifgs"
        mkdir -p $ifgdir.backup
        mv $ifgdir/* $ifgdir.backup/.
    fi
};


function correct_ifg_tides_public() {
  if [ "$#" == "3" ]; then
    local frame=$1
    local ifg=$2
    local ext=$3
    track=`echo $frame | cut -c -3 | sed 's/^0//' | sed 's/^0//'`
    ifgpath=$LiCSAR_public/$track/$frame/interferograms/$ifg
    infile=$ifgpath/$ifg.geo.$ext.tif
    if [ ! -f $infile ]; then
      echo "interferogram does not exist"
      return 0
    fi
    outfile=$ifgpath/$ifg.geo.$ext.notides.tif
    if [ -f $outfile ]; then
      echo "tide corrected ifg "$ifg" already exists"
      return 0
    fi
    date1=`echo $ifg | cut -d '_' -f1`
    date2=`echo $ifg | cut -d '_' -f2`
    create_LOS_tide $frame $date1
    create_LOS_tide $frame $date2
    tided1=$LiCSAR_public/$track/$frame/epochs/$date1/$date1'.tide.geo.tif'
    tided2=$LiCSAR_public/$track/$frame/epochs/$date2/$date2'.tide.geo.tif'
    #wavelength of Sentinel-1
    #wavelength=0.055465763
    #mcoef=-4*PI/$wavelength
    #mcoef=-226.56
    echo "correcting the interferogram to "$outfile
    gmt grdmath $tided2 $tided1 SUB -226.56 MUL = $outfile.onlytide.nc
    if [ $ext == "unw" ]; then
     gmt grdmath -N $infile'=gd:Gtiff+n0' 0 NAN $tided2 $tided1 SUB -226.56 MUL SUB = $outfile.nc #'=gd:Gtiff'
     gmt grdmath $outfile.nc $outfile.nc MEDIAN SUB = $outfile=gd:Gtiff
     #demedian_unw.py $outfile'=gd:Gtiff'
     #create_preview_unwrapped $outfile
    else
     gmt grdmath -N $infile'=gd:Gtiff+n0' 0 NAN $tided2 $tided1 SUB -226.56 MUL SUB WRAP = $outfile'=gd:Gtiff'
     create_preview_wrapped $outfile
    fi
  else
    echo "Usage: correct_ifg_tides_public frame_id ifg ext"
    echo "(ifg must be yyyymmdd_yyyymmdd)"
    echo "(ext should be either diff_pha or diff_unfiltered_pha)"
    return 0
  fi
};


function get_burstsno_frame(){
 if [ -z $1 ]; then echo "just add frame"; return 0; fi
 frame=$1
 bursts=`echo $frame | cut -d '_' -f3`
 burstsnum=`echo ${bursts:0:2}+${bursts:2:2}+${bursts:4:2} | bc`
 echo $burstsnum
}

function get_master(){
 tdir=`pwd`
 if [ ! -z $1 ]; then
  cdproc $1
 fi
 #this should be run from a frame folder
 if [ ! -d geo ]; then
  #echo "please use get_master function only in frame folder, with geo"
  cd $tdir
  return 0
 else
  m=`ls geo/20??????.hgt 2>/dev/null| cut -d '/' -f2 | cut -d '.' -f1`
  if [ -z $m ]; then
   #echo "error - missing hgt file in geo folder"
   cd $tdir
   return 0
  else
   echo $m
  fi
 fi
 cd $tdir
};


function init_LiCSAR(){
  # RAW_DIR variable has to be defined in the launch script!!
  while getopts ":m:p:" OPTION ; do
    case $OPTION in
    m)    master="$OPTARG";;
    p)    polygon_file="$OPTARG";;
    esac
  done
  shift $((OPTIND-1))
  if [[ "$master" == "-m" || "$master" == "-p" || -z "${master}" ]]; then master=0; fi
  if [[ "$master" == "-p" || -z "${polygon_file}" ]]; then polygon_file=0; fi
  
  # If master and polygon exists proceed with the cleaning of zipfile list
  if [[ "${master}" != "0" && "${polygon_file}" != "0" ]]; then
    echo "Master date and Polygon file available... Computing which zipfiles are in bounds of polygon"
    # generate the list of zip files to be linked [Future versions should grab this info from the database]
    createzipfileslist ${RAW_DIR} $polygon_file 
    if [ -z "$burstid_file" ]; then
        burstid_file=$(basename $polygon_file .xy)_burst_ids.txt  
        master2burstlist ${RAW_DIR} ${master} ${polygon_file} ${burstid_file} # Generate the burst centers list based on the polygon_file       
      fi
  else
    if [[ "${master}" == "0" && "${polygon_file}" != "0" ]]; then # If master info is ok but polygon is missing
      echo "Polygon seems to be ok, compute master based on polygon: $polygon_file master: $master "
      echo "Computing which zipfiles are in bounds of polygon..."
      # generate the list of zip files to be linked [Future versions should grab this info from the database]
      echo "createzipfileslist ${RAW_DIR} $polygon_file "
      createzipfileslist ${RAW_DIR} $polygon_file 
      echo "createzipfileslist ${RAW_DIR} $polygon_file "
      # Check if master date exists, if not assign one. That makes way easier case handling internally in LiCSAR
      check_master zipfile_inbounds.list $master ; master=`cat master_date.txt`  
      if [ -z "$burstid_file" ]; then
        burstid_file=$(basename $polygon_file .xy)_burst_ids.txt  
        master2burstlist ${RAW_DIR} ${master} ${polygon_file} ${burstid_file} # Generate the burst centers list based on the polygon_file       
      fi
    elif [[  "${master}" != "0" && "${polygon_file}" == "0" ]]; then # If polygon is ok but master info is missing
      echo "Master seems to be ok, compute polygon based on master: $master polygon: $polygon_file"
      polygon_file=auto_polygon.xy
      burstid_file=$(basename $polygon_file .xy)_burst_ids.txt      
      echo "Creating polygon (${polygon_file}) and burstid files (${burstid_file})"
      boundlims2polygon ${RAW_DIR} ${master} ${polygon_file} # Generate the list of coordinates of the limits we want to include
      master2burstlist ${RAW_DIR} ${master} ${polygon_file} ${burstid_file} # Generate the burst centers list based on the polygon_file
      createzipfileslist ${RAW_DIR} $polygon_file 
    else  # Else there is no info on the master and polygon, then abort
      echo "Either master or polygon_file do not exist. We have to select / generate these information"
      exit 0;
    fi
  fi
  echo "Master date: $master ; and Polygon file: $polygon_file "
}; export -f init_LiCSAR


function boundlims2polygon(){
  #######################################################
  # create a file with a list of coordinates based on some bounding limits
  #######################################################
  if [ "$#" == "5" ]; then
    local west=$1
    local east=$2
    local south=$3
    local north=$4  
    local outfile=$5
    rm -f $outfile
    echo "$west $north" >> $outfile # upper left
    echo "$east $north" >> $outfile # upper right
    echo "$east $south" >> $outfile # lower right
    echo "$west $south" >> $outfile # lower left
    echo "$west $north" >> $outfile # upper left
  elif [ "$#" == "3" ]; then
    local pathmaster=$1
    local master=$2
    local outfile=$3
    for i in `ls ${pathmaster}/*${master}*.zip`; do
      #echo "zipSAFE_2_bounding_box.py $i >> coordinates.$$" 
      zipSAFE_2_SLCbb.py $i >> coordinates.$$
    done
    SLCbb_minmax.py coordinates.$$ > coordinates
    local west=`awk '{print $1}' coordinates`; 
    local east=`awk '{print $2}' coordinates`;
    local south=`awk '{print $3}' coordinates`;
    local north=`awk '{print $4}' coordinates`;
    echo "$west $north" > $outfile # upper left
    echo "$east $north" >> $outfile # upper right
    echo "$east $south" >> $outfile # lower right
    echo "$west $south" >> $outfile # lower left
    echo "$west $north" >> $outfile # upper left    
    rm coordinates*
  else
    echo "We cannot find the information needed to generate a polygon file"
    exit 0
  fi
}; export -f boundlims2polygon # End function definition

function master2burstlist(){
  #######################################################
  # create a file with a list of coordinates based on some bounding limits
  #######################################################
  if [ "$#" == "4" ]; then
    local pathmaster=$1
    local master=$2
    local polygon_file=$3
    local outfile=$4
    for i in `ls ${pathmaster}/*${master}*.zip`; do
      #echo "zipSAFE_2_BurstList.py $i" 
      zipSAFE_2_BurstsList.py $i #>> coordinates.$$
    done
    for i in `seq 1 3`; do
      cat *T*_IW${i}.burstlist > ${master}_IW${i}.burstlist
      BurstListCoordInPolyFile.py ${master}_IW${i}.burstlist $polygon_file > IW${i}.burstlist
    done
    echo "IW1 `wc -lc IW1.burstlist | awk '{print $1}'` IW2 `wc -lc IW2.burstlist | awk '{print $1}'` IW3 `wc -lc IW3.burstlist | awk '{print $1}'`" > $outfile
    awk '{print "burstID", $1, $2}' IW1.burstlist >> $outfile
    awk '{print "burstID", $1, $2}' IW2.burstlist >> $outfile
    awk '{print "burstID", $1, $2}' IW3.burstlist >> $outfile
    #rm -f *IW*.burstlist 
  else
    echo "We cannot find the information needed to generate burst list for the master date image"
    exit 0
  fi
}; export -f master2burstlist # End function definition


function createzipfileslist(){
  RAW_DIR=$1
  polygon_file=$2
  # RAW_DIR is an environmental variable created in LiCSAR_config_launch
  # generate the list of zip files to be linked [Future versions should grab this info from the database]
  echo "ls ${RAW_DIR}/S1?_*.zip > zipfile.list"
  ls ${RAW_DIR}/S1?_*.zip > zipfile.list
  # Check if each file in the ziplist falls within the bounding limits of the frame
  check_ziplist_inbounds.sh zipfile.list $polygon_file > zipfile_inbounds.list
}; export -f createzipfileslist

function sortarray(){
# Function helps to sort a list of SLC zip files according to increasing dates
# This is motivated that dates for dual and single polarization are differently sorted
  declare -a inarray=("${!1}")
  slcdir=$2
  rm -f unsortedarray.list
  nslcs=${#inarray[@]}
  maxi=`echo $nslcs | awk '{print $1-1}'`
  for i in `seq 0 ${maxi}`; do
    basename ${inarray[${i}]} >> unsortedarray.list ; 
  done ; 
  # Sort list based on positions from 18 to 25
  sort -k1.18,1.25 unsortedarray.list > sortedarray.list
  awk '{print "'${slcdir}'"$1}' sortedarray.list > sortedarray2.list
  slclist=( `cat sortedarray2.list` )
  rm -f sortedarray.list sortedarray2.list unsortedarray.list
}; export -f sortarray # End function definition

function sortarray_slc_yyyymmdd(){
# Function helps to sort a list of SLC zip files according to increasing dates
# This is motivated that dates for dual and single polarization are differently sorted
  declare -a inarray=("${!1}")
  slcdir=$2
  rm -f unsortedarray.list
  nslcs=${#inarray[@]}
  maxi=`echo $nslcs | awk '{print $1-1}'`
  for i in `seq 0 ${maxi}`; do
    basename ${inarray[${i}]} >> unsortedarray.list ; 
  done ; 
  # Sort list based on positions from 18 to 25
  sort -k1.18,1.25 unsortedarray.list > sortedarray.list
  awk '{print "'${slcdir}'"$1}' sortedarray.list > sortedarray2.list
  slclist=( `cat sortedarray2.list` )
  rm -f sortedarray.list sortedarray2.list unsortedarray.list
}; export -f sortarray_slc_yyyymmdd # End function definition

function cp_OrbitFile(){
# I think it is fixed to find:
# Restituted orbits need to specify the time to the hour (time of data - 3 hours)
# Precise orbits need to specify the time to the day (time of data - 1 day)

  local imageID=$1
  local srcOrbdir=$2
  local mvdir=$3
  yyyymmdd=`basename ${imageID} | awk '{print substr($1,18,8)}'`
  ABCD=`basename ${imageID} | awk '{print substr($1,3,1)}'`
  # For Precise orbits we need '$yyyymmdd'-1
  yyyymmddmin1d=`date "--date=${yyyymmdd} -1 day" +%Y%m%d` 
  hh=`basename ${imageID} | awk '{print substr($1,27,2)}'`
  mm=`basename ${imageID} | awk '{print substr($1,29,2)}'`
  ss=`basename ${imageID} | awk '{print substr($1,31,2)}'`

  # For Restituted orbits we need to be within a certain time window (~3 hours)
  #echo "--date=${yyyymmdd} ${hh}:${mm}:${ss}"
  yyyymmdd2sec=`date "--date=${yyyymmdd} ${hh}:${mm}:${ss}" +%s` # Convert to seconds since 1970-01-01 00:00:00 UTC
  
  # Find the Precise or Restituted orbit file
  if [ -f "`ls $srcOrbdir/POEORB/S1${ABCD}_*V${yyyymmddmin1d}*.EOF`" ]; then
    orbit_file=`ls $srcOrbdir/POEORB/S1${ABCD}_*V${yyyymmddmin1d}*.EOF`
    cp $orbit_file $mvdir
    echo "    SLC: $yyyymmdd.slc has a PRECISE orbit file: $orbit_file"
  else 
    # Loop through them and find the matching one
    for fileID in `ls $srcOrbdir/RESORB/S1${ABCD}_*V${yyyymmdd}*.EOF`; do
      hhdate1=`basename $fileID | awk '{print substr($1,52,2)}'`
      mmdate1=`basename $fileID | awk '{print substr($1,54,2)}'`
      ssdate1=`basename $fileID | awk '{print substr($1,56,2)}'`
      hhdate2=`basename $fileID | awk '{print substr($1,68,2)}'`
      mmdate2=`basename $fileID | awk '{print substr($1,70,2)}'`
      ssdate2=`basename $fileID | awk '{print substr($1,72,2)}'`
      resorbbegin=`date "--date=${yyyymmdd} ${hhdate1}:${mmdate1}:${ssdate1}" +%s`      
      resorbend=`date "--date=${yyyymmdd} ${hhdate2}:${mmdate2}:${ssdate2}" +%s`      
      if [[ ( "$resorbbegin" -lt "$yyyymmdd2sec" ) && ( "$resorbend" -gt "$yyyymmdd2sec" ) ]]; then
        orbit_file=$fileID
        cp $orbit_file $mvdir
        echo "    SLC: $yyyymmdd.slc has a RESTITUTED orbit file: $orbit_file"
        break # Once we find one orbit file we skip it
      fi
    done
  fi
  if [ -z $orbit_file ]; then
    echo "SLC: $yyyymmdd.slc has NO ORBIT FILE. Now, we try to download it"
    fetchOrbit_python27.py -i ${imageID} -o SLC/$yyyymmdd/ 
  fi
}; export -f cp_OrbitFile # End function definition

########################################################################
########################################################################
# Apply the orbits to the slc files
function UpdateOrbits(){
  local yyyymmdd=$1
  local ORBs_DIR=$2
  
  # Download precise orbits information (available from https://qc.sentinel1.eo.esa.int/aux_poeorb/) and update .slc.par files
  # For each product, let's download the orbit files (sometimes it does more than one, if multiple slices from same date are in SLC/ folder)
  for i in `ls SLC/${yyyymmdd}/S1?*.zip`; do
    zipfilename=`basename ${i}`
    cp_OrbitFile ${zipfilename} $ORBs_DIR SLC/${yyyymmdd}
  done
  
  # For each merged product applies the orbit file
  # Find the Precise or Restituted orbit file
  yyyymmddmin=`date "--date=${yyyymmdd} -1 day" +%Y%m%d`
  if [ -f "`ls SLC/${yyyymmdd}/S1?_*V${yyyymmdd}*.EOF`" ]; then
    orbit_file=`ls SLC/${yyyymmdd}/S1?_*V${yyyymmdd}*.EOF`
    echo "    SLC: $yyyymmdd.slc has a RESTITUTED orbit file: $orbit_file"
    for img in `ls SLC/${yyyymmdd}/${yyyymmdd}*.sl*.par` ; do
      echo "      S1_OPOD_vec ${img} $orbit_file "
      S1_OPOD_vec ${img} $orbit_file 
    done
  elif [ -f "`ls SLC/${yyyymmdd}/S1?_*V${yyyymmddmin}*.EOF`" ]; then
    orbit_file=`ls SLC/${yyyymmdd}/S1?_*V${yyyymmddmin}*.EOF`
    echo "    SLC: $yyyymmdd.slc has a PRECISE orbit file: $orbit_file"
    for img in `ls SLC/${yyyymmdd}/${yyyymmdd}*.sl*.par` ; do
      echo "      S1_OPOD_vec ${img} $orbit_file "      
      S1_OPOD_vec ${img} $orbit_file 
    done
  else
    echo "SLC: $yyyymmdd.slc has NO ORBIT FILE"
  fi
}; export -f UpdateOrbits  # End function definition

########################################################################
########################################################################
# Pad with zeros the input slc
function pad_SLCs(){
  # Image date to be processed if padding is necessary
  local yyyymmdd=$1
  # Input variables IW1 and IW2 indicates which subswaths to be used 
  # IW1-to-IW2, IW2-to-IW3 or IW1-to-IW3, if not provided assumes IW1-to-IW3 
  if [ -z "$2" ]; then local IW1=1; else local IW1=$2; fi 
  if [ -z "$3" ]; then local IW2=3; else local IW2=$3; fi   
  ##rjw ## Before merging, pad bursts with extra lines if necessary ##
  echo "  Checking num. burst lines tally for all ${yyyymmdd} slc files"
  for n in `seq ${IW1} ${IW2}`; do
    maxburst_smp=`grep lines_per_burst SLC/${yyyymmdd}/${yyyymmdd}T*IW$n.slc.TOPS_par | awk 'n < $2 {n=$2}END{print n}'`
    for yyyymmddhhmmss in `ls -d SLC/${yyyymmdd}/${yyyymmdd}T*IW$n.slc | awk '{print substr($1,14,15)}'`; do
      rng_smp=`grep range_samples SLC/${yyyymmdd}/$yyyymmddhhmmss*IW$n.slc.par | awk '{print $2}'`
      azm_smp=`grep azimuth_lines SLC/${yyyymmdd}/$yyyymmddhhmmss*IW$n.slc.par | awk '{print $2}'`
      nbursts=`grep number_of_bursts SLC/${yyyymmdd}/$yyyymmddhhmmss*IW$n.slc.TOPS_par | awk '{print $2}'`
      burst_smp=`grep lines_per_burst SLC/${yyyymmdd}/$yyyymmddhhmmss*IW$n.slc.TOPS_par | awk '{print $2}'`
      diff_smp=`echo $maxburst_smp $burst_smp | awk '{print $1 - $2}'`
      if [ $diff_smp -ne "0" ]; then
        echo "SLC" SLC/$yyyymmdd/$yyyymmddhhmmss*IW$n.slc "has" $diff_smp "less lines per burst than SLC with maximum number, making correction..."
        burst_bytes=`echo $burst_smp $rng_smp | awk '{print $1*$2*4}'`
        pad_bytes=`echo $diff_smp $rng_smp | awk '{print $1*$2*4}'` 
        # assumes short (4bytes per pixel)

        echo "Splitting file" SLC/$yyyymmdd/$yyyymmddhhmmss*IW$n.slc "into" $nbursts "bursts of" $burst_bytes "bytes each..."
        split -b $burst_bytes -d SLC/$yyyymmdd/$yyyymmddhhmmss*IW$n.slc SLC/$yyyymmdd/burst
        tail -c $pad_bytes SLC/$yyyymmdd/burst00 > SLC/$yyyymmdd/pad01

        echo "Stitching new SLC with padding between bursts"
        catlist=`ls SLC/$yyyymmdd/burst* | awk '{ORS=" "; print $1, "SLC/'"$yyyymmdd"'/pad01"}'`
        cat $catlist > SLC/$yyyymmdd/$yyyymmddhhmmss*IW$n.pad.slc 
        rm -f SLC/$yyyymmdd/burst* SLC/$yyyymmdd/pad01
        mv SLC/$yyyymmdd/$yyyymmddhhmmss*IW$n.pad.slc SLC/$yyyymmdd/$yyyymmddhhmmss*IW$n.slc

        echo "Editing" SLC/$yyyymmdd/$yyyymmddhhmmss*IW$n.slc.par "and" SLC/$yyyymmdd/$yyyymmddhhmmss*IW$n.slc.TOPS_par
        newazm_smp=`echo $azm_smp $nbursts | awk '{print $1*$2}'`
        sed -i 's/azimuth_lines: *'$azm_smp'/azimuth_lines:                 '$newazm_smp'/g' SLC/$yyyymmdd/$yyyymmddhhmmss*IW$n.slc.par 
        sed -i 's/lines_per_burst: *'$burst_smp'/lines_per_burst:   '$maxburst_smp'/g'  SLC/$yyyymmdd/$yyyymmddhhmmss*IW$n.slc.TOPS_par 
      fi
    done
  done
  ############# rjw
}; export -f pad_SLCs # End function definition

########################################################################
########################################################################
# 
function createSLCtab(){
  local inimage=$1
  local ending=$2
  if [ -z "$3" ]; then local IW1=1; else local IW1=$3; fi 
  if [ -z "$4" ]; then local IW2=3; else local IW2=$4; fi   
  for i in `seq ${IW1} ${IW2}`; do
    echo "${inimage}.IW${i}.${ending} ${inimage}.IW${i}.${ending}.par ${inimage}.IW${i}.${ending}.TOPS_par"
  done
}; export -f createSLCtab # End function definition

########################################################################
########################################################################
# 
function merge_S1cat(){
  local slcdir=$1
  local date1=$2
  local date2=$3
  local date=$4
  if [ -z "$5" ]; then local IW1=1; else local IW1=$5; fi 
  if [ -z "$6" ]; then local IW2=3; else local IW2=$6; fi  
  createSLCtab $slcdir/$date1 slc ${IW1} ${IW2} > SLC_tab1 
  createSLCtab $slcdir/$date2 slc ${IW1} ${IW2} > SLC_tab2
  createSLCtab $slcdir/$date  slc ${IW1} ${IW2} > SLC_tab3
  SLC_cat_S1_TOPS SLC_tab1 SLC_tab2 SLC_tab3
}; export -f merge_S1cat # End function definition

########################################################################
########################################################################
# Move files from a long name to a shorter name
function mvSLC_yyyymmddhhmmss2yyyymmdd(){
  local yyyymmddhhmmss=$1
  local yyyymmdd=$2
  if [ -z "$3" ]; then local IW1=1; else local IW1=$3; fi 
  if [ -z "$4" ]; then local IW2=3; else local IW2=$4; fi   
  for i in `seq ${IW1} ${IW2}`; do
    mv SLC/${yyyymmdd}/${yyyymmddhhmmss}.IW${i}.slc           SLC/${yyyymmdd}/${yyyymmdd}.IW${i}.slc 
    mv SLC/${yyyymmdd}/${yyyymmddhhmmss}.IW${i}.slc.par       SLC/${yyyymmdd}/${yyyymmdd}.IW${i}.slc.par 
    mv SLC/${yyyymmdd}/${yyyymmddhhmmss}.IW${i}.slc.TOPS_par  SLC/${yyyymmdd}/${yyyymmdd}.IW${i}.slc.TOPS_par
  done
}; export -f mvSLC_yyyymmddhhmmss2yyyymmdd # End function definition

########################################################################
########################################################################
# Move files from a cropped version to normal version (Merging step)
function mvSLCcrop2SLC(){
  local iput=$1
  local oput=$2
  if [ -z "$3" ]; then local IW1=1; else local IW1=$3; fi 
  if [ -z "$4" ]; then local IW2=3; else local IW2=$4; fi   
  for i in `seq ${IW1} ${IW2}`; do
    mv ${iput}.IW${i}.slc           ${oput}.IW${i}.slc 
    mv ${iput}.IW${i}.slc.par       ${oput}.IW${i}.slc.par 
    mv ${iput}.IW${i}.slc.TOPS_par  ${oput}.IW${i}.slc.TOPS_par
  done
}; export -f mvSLCcrop2SLC # End function definition

########################################################################
########################################################################
# Crop files according to list of common bursts
function mk_crop(){
  local image1=$1
  if [ -z "$2" ]; then local IW1=1; else local IW1=$2; fi 
  if [ -z "$3" ]; then local IW2=3; else local IW2=$3; fi   
  for i in `seq ${IW1} ${IW2}`; do
    MinBurstIW=`awk '{print $1}' ${image1}_IW${i}.commonburst | sort -n | head -1`;
    MaxBurstIW=`awk '{print $2}' ${image1}_IW${i}.commonburst | sort -n | tail -1`;  
    createSLCtab ${image1} slc ${i} ${i} > SLCin_tab
    createSLCtab ${image1}.crop slc ${i} ${i} > SLCincrop_tab
    echo "  Cropping image $image1 IW${i}: SLC_copy_S1_TOPS SLCin_tab SLCincrop_tab 1 $MinBurstIW 1 $MaxBurstIW"
    SLC_copy_S1_TOPS SLCin_tab SLCincrop_tab 1 $MinBurstIW 1 $MaxBurstIW 
  done
  # Mosaic/Merge cropped subswaths
  createSLCtab ${image1}.crop slc ${IW1} ${IW2} > SLCincropped_tab
  echo "  SLC_mosaic_S1_TOPS SLCincropped_tab ${image1}.slc ${image1}.slc.par $rlks $azlks " 
  SLC_mosaic_S1_TOPS SLCincropped_tab ${image1}.slc ${image1}.slc.par $rlks $azlks  
  mvSLCcrop2SLC ${image1}.crop ${image1} ${IW1} ${IW2} 

  # Update the ${master}_IW${i}.burstlist files
  for i in `seq ${IW1} ${IW2}`; do
    MinBurstIW=`awk '{print $1}' ${image1}_IW${i}.commonburst | sort -n | head -1`;
    MaxBurstIW=`awk '{print $2}' ${image1}_IW${i}.commonburst | sort -n | tail -1`;  
    #echo $i $MinBurstIW $MaxBurstIW ${image1}_IW${i}.burstlist
    awk 'NR>='$MinBurstIW'&&NR<='$MaxBurstIW' {print $1,$2,NR-'$MinBurstIW'+1}' ${image1}_IW${i}.burstlist > temp
    mv temp ${image1}_IW${i}.burstlist
  done
  
  # Cleaning
  rm -f SLCin_tab SLCincrop_tab SLCincropped_tab 
}; export -f mk_crop # End function definition

########################################################################
########################################################################
# Crop files according to list of common bursts
function mk_recrop(){
  local image1=$1
  local outimage1=$2
  if [ -z "$3" ]; then local IW1=1; else local IW1=$3; fi 
  if [ -z "$4" ]; then local IW2=3; else local IW2=$4; fi   
  for i in `seq ${IW1} ${IW2}`; do
    MinBurstIW=`awk '{print $1}' ${outimage1}_IW${i}.commonburst | sort -n | head -1`;
    MaxBurstIW=`awk '{print $2}' ${outimage1}_IW${i}.commonburst | sort -n | tail -1`;  
    createSLCtab ${image1} slc ${i} ${i} > SLCin_tab
    createSLCtab ${outimage1}.crop slc ${i} ${i} > SLCincrop_tab
    echo "  Cropping image $image1 IW${i}: SLC_copy_S1_TOPS SLCin_tab SLCincrop_tab 1 $MinBurstIW 1 $MaxBurstIW"
    SLC_copy_S1_TOPS SLCin_tab SLCincrop_tab 1 $MinBurstIW 1 $MaxBurstIW 
  done
  # Mosaic/Merge cropped subswaths
  createSLCtab ${outimage1}.crop slc ${IW1} ${IW2} > SLCincropped_tab
  echo "  SLC_mosaic_S1_TOPS SLCincropped_tab ${outimage1}.slc ${outimage1}.slc.par $rlks $azlks " 
  SLC_mosaic_S1_TOPS SLCincropped_tab ${outimage1}.slc ${outimage1}.slc.par $rlks $azlks  
  mvSLCcrop2SLC ${outimage1}.crop ${outimage1} ${IW1} ${IW2} 

  # Cleaning
  rm -f SLCin_tab SLCincrop_tab SLCincropped_tab 
}; export -f mk_recrop # End function definition


function check_missing_bursts(){ 
  local polygon_file=$1 
  local burstid_file=$2 
  local locdir=$3 
  local yyyymmdd=$4
  if [ -z "$5" ]; then local IW1=1; else local IW1=$5; fi 
  if [ -z "$6" ]; then local IW3=3; else local IW3=$6; fi 
  # Identify which bursts in master are inside the polygon, before starting to merge
  rm -f ${locdir}/*_zipfiles.txt ${locdir}/*_zipfiles_short.txt
    
  # Check the zipfiles 
  burstsinframe2.sh ${polygon_file} ${burstid_file} ${locdir} ${yyyymmdd} 
  nb_inframe=`wc -l < $burstid_file | awk '{print $1-1}'`
  nb_inimage=`cat SLC/${yyyymmdd}/${yyyymmdd}_zipfiles_short.txt | grep -o NA | wc -l`
  
  # Check if the number of bursts missing is larger than zero
  if [ "$nb_inimage" == "0" ]; then 
    echo $nb_inimage
  else
    # First loop from IW1 to IW3 
    for i in `seq $IW1 $IW3`; do
      awk '{if ($4=="'${i}'") print $0}' SLC/${yyyymmdd}/${yyyymmdd}_zipfiles_short.txt > tmp1
      # Convert fifth columns in 1 or 0 (awk '{if ($5=="NA") print 0; else print 1}' tmp1)
      # Then remove repeated consecutive lines (awk '{sub(".*"$1,$1)}1' | uniq)
      # Finally transpose the results to find/match a pattern
      awk '{if ($5=="NA") print 0; else print 1}' tmp1 | awk '{sub(".*"$1,$1)}1' | uniq | paste -sd" " | awk '{ if (/1 0 1/) { print 1} else {print 0} }' >> tmpvals
    done
    datagaps=`awk '{s+=$1} END {print s}' tmpvals`
    if [ "$datagaps" == "0" ]; then 
      echo $nb_inimage
    else
      echo "-1"
    fi
    rm -f tmpvals tmp1
  fi  
}; export -f check_missing_bursts

########################################################################
########################################################################
# Uncompress zip files
function uncompressSAFE(){
  local slcdate=$1
  local slcdir=SLC/${slcdate}
  cd ${slcdir}
  for i in `ls *.zip`; do 
    echo "  Uncompressing for date: $slcdate file: $i"
    #jar -xf $i ; # It is necessary to move where the zip softlinks are to uncompress them
    7za x $i ; # It is necessary to move where the zip softlinks are to uncompress them
  done
  cd ../..
}; export -f uncompressSAFE

########################################################################
########################################################################
# Read a geotiff file into GAMMA format
function readSAFE(){
  local slcdate=$1
  local slcdir=SLC/${slcdate}
  for productID in `ls -d ${slcdir}/S1?*.SAFE/ ` ; do
    echo "  Reformating $productID file"
    yyyymmddhhmmss=`echo ${productID} | awk '{print substr($1,31,15)}'`
    yyyymmdd=${yyyymmddhhmmss:0:8}
    geotiff_vv_iw1=${productID}/*/s1?-iw1-slc-vv-${yyyymmdd}*.tiff
    geotiff_vv_iw2=${productID}/*/s1?-iw2-slc-vv-${yyyymmdd}*.tiff
    geotiff_vv_iw3=${productID}/*/s1?-iw3-slc-vv-${yyyymmdd}*.tiff
    ann_vv_iw1=${productID}/*/s1?-iw1-slc-vv-${yyyymmdd}*.xml
    ann_vv_iw2=${productID}/*/s1?-iw2-slc-vv-${yyyymmdd}*.xml
    ann_vv_iw3=${productID}/*/s1?-iw3-slc-vv-${yyyymmdd}*.xml
    cal_vv_iw1=${productID}/*/*/calibration-s1?-iw1-slc-vv-${yyyymmdd}*.xml
    cal_vv_iw2=${productID}/*/*/calibration-s1?-iw2-slc-vv-${yyyymmdd}*.xml
    cal_vv_iw3=${productID}/*/*/calibration-s1?-iw3-slc-vv-${yyyymmdd}*.xml
    noise_vv_iw1=${productID}/*/*/noise-s1?-iw1-slc-vv-${yyyymmdd}*.xml
    noise_vv_iw2=${productID}/*/*/noise-s1?-iw2-slc-vv-${yyyymmdd}*.xml
    noise_vv_iw3=${productID}/*/*/noise-s1?-iw3-slc-vv-${yyyymmdd}*.xml
   
    # Check again for short integer rather than float32 complex for the .slc files
    par_S1_SLC ${geotiff_vv_iw1} ${ann_vv_iw1} ${cal_vv_iw1} ${noise_vv_iw1} ${slcdir}/${yyyymmddhhmmss}.IW1.slc.par ${slcdir}/${yyyymmddhhmmss}.IW1.slc ${slcdir}/${yyyymmddhhmmss}.IW1.slc.TOPS_par 1 60 
    par_S1_SLC ${geotiff_vv_iw2} ${ann_vv_iw2} ${cal_vv_iw2} ${noise_vv_iw2} ${slcdir}/${yyyymmddhhmmss}.IW2.slc.par ${slcdir}/${yyyymmddhhmmss}.IW2.slc ${slcdir}/${yyyymmddhhmmss}.IW2.slc.TOPS_par 1 60 
    par_S1_SLC ${geotiff_vv_iw3} ${ann_vv_iw3} ${cal_vv_iw3} ${noise_vv_iw3} ${slcdir}/${yyyymmddhhmmss}.IW3.slc.par ${slcdir}/${yyyymmddhhmmss}.IW3.slc ${slcdir}/${yyyymmddhhmmss}.IW3.slc.TOPS_par 1 60 
    
    # Delete folder .SAFE to save disk storage
    rm -rf $productID
  done
}; export -f readSAFE

########################################################################
########################################################################
# Merge all available SLCs files in a SLC/yyyymmdd directory to generate a single SLC TOPS file
function mergeSLC(){
  local slcdate=$1
  local rlks=$2 
  local azlks=$3
  local slcdir=SLC/${slcdate}
  # Merge the burst list files
  echo "  Merging files for date: $slcdate [if only one file skip merging, only rename it] "
  rm -f ${slcdir}/${slcdate}_IW*.burstlist
  cat ${slcdir}/${slcdate}T*_IW1.burstlist > ${slcdir}/${slcdate}_IW1.burstlist
  cat ${slcdir}/${slcdate}T*_IW2.burstlist > ${slcdir}/${slcdate}_IW2.burstlist
  cat ${slcdir}/${slcdate}T*_IW3.burstlist > ${slcdir}/${slcdate}_IW3.burstlist
  # Check for number of SLCs per date
  listslc=( `ls ${slcdir}/${slcdate}T??????.IW1.slc | awk '{print substr($1,14,15)}'` );
  nSLCs=${#listslc[@]};
  if [ "$nSLCs" == "1" ]; then
    echo "    Only 1 SLC per date is available [rename] "
    yyyymmddhhmmss=`ls ${slcdir}/${slcdate}T??????.IW1.slc | awk '{print substr($1,14,15)}'`  
    mvSLC_yyyymmddhhmmss2yyyymmdd ${yyyymmddhhmmss} ${slcdate} 
    createSLCtab ${slcdir}/${slcdate} slc ${IW1} ${IW3} > SLC_tab
    # Create a mosaic SLC from the three subswath bursted SLCs
    SLC_mosaic_S1_TOPS SLC_tab ${slcdir}/${slcdate}.slc ${slcdir}/${slcdate}.slc.par $rlks $azlks 
    rm -f SLC_tab
  else
    echo "    ${nSLCs} SLCs in this date are available for merging"
    # Loop from first element in array to last element minus 2 to merge in pairs
    nslcmintwo=`echo $nSLCs | awk '{print $1-2}'`;
    for i in `seq 0 $nslcmintwo`; do # Indexes in bash arrays starts with 0!!!
      echo ${listslc[*]}
      j=`echo $i | awk '{print $1+1}'`
      SLC1=${listslc[${i}]}; # First image to concatenate
      SLC2=${listslc[${j}]}; # Second image to concatenate       
      if [ "$i" -lt "$nslcmintwo" ]; then 
        echo "    Inputs are: $SLC1 $SLC2     Output: ${slcdate}.${i}.${j} "
        listslc[${j}]=${slcdate}.${i}.${j} # Substitute the jth element in array for the new name
        # It merges and rename products as yyyymmdd format
        merge_S1cat $slcdir $SLC1 $SLC2 ${slcdate}.${i}.${j} # Concatenate pairs and generates SLC tables (SLC_tabs) 
        rm -f SLC_tab1 SLC_tab2 SLC_tab3         
      else
        echo "    Inputs are: $SLC1 $SLC2     Output: ${slcdate} "
        # It merges and rename products as yyyymmdd format
        merge_S1cat $slcdir $SLC1 $SLC2 $slcdate # Concatenate pairs and generates SLC tables (SLC_tabs) 
        SLC_mosaic_S1_TOPS SLC_tab3 ${slcdir}/${slcdate}.slc ${slcdir}/${slcdate}.slc.par $rlks $azlks 
        rm -f SLC_tab1 SLC_tab2 SLC_tab3         
      fi 
    done # End of the sequential concatenatation merging
  fi
}; export -f mergeSLC

########################################################################
########################################################################
# Generate a multilooked version of the mosaic in a TOPS file (optionally a raster image)
function multilookSLC(){
  local slcdate=$1
  local rlks=$2 
  local azlks=$3
  local plotme=$4
  local slcdir=$5
  if [ -z $slcdir ]; then slcdir=SLC/$slcdate; fi
  # Generate a quicklook to inspect everything went correctly
  echo "    Multilooking image ${slcdate} multilook factor $rlks [range] $azlks [azimuth] "
  multi_look $slcdir/${slcdate}.slc $slcdir/${slcdate}.slc.par $slcdir/${slcdate}.slc.mli $slcdir/${slcdate}.slc.mli.par $rlks $azlks > /dev/null 2>&1 #>> $logfile 
  # Create a thumbnail image to inspect the mosaic      
  if [ "$plotme" == "1" ]; then 
    wid=`grep range_samples: $slcdir/${slcdate}.slc.mli.par | awk '{print $2}'` 
    reducfac=`echo $wid | awk '{if(int($1/1000) > 1) print int($1/1000); else print 1}'` 
    raspwr $slcdir/${slcdate}.slc.mli $wid - - $reducfac $reducfac - - - $slcdir/${slcdate}.slc.mli.bmp 0 > /dev/null 2>&1 # >> $logfile 
  fi
}; export -f multilookSLC

########################################################################
########################################################################
# Generate a multilooked version of the mosaic in a Resampled TOPS file (optionally a raster image)
function multilookRSLC(){
  local slcdate=$1
  local rlks=$2 
  local azlks=$3
  local plotme=$4
  local rslcdir=$5
  if [ -z $rslcdir ]; then rslcdir=RSLC/$slcdate; fi
  # Generate a quicklook to inspect everything went correctly
  echo "    Multilooking image ${slcdate} multilook factor $rlks [range] $azlks [azimuth] "
  multi_look ${rslcdir}/${slcdate}.rslc ${rslcdir}/${slcdate}.rslc.par ${rslcdir}/${slcdate}.rslc.mli ${rslcdir}/${slcdate}.rslc.mli.par $rlks $azlks > /dev/null 2>&1 #>> $logfile 
  # Create a thumbnail image to inspect the mosaic      
  if [ "$plotme" == "1" ]; then 
    wid=`grep range_samples: ${rslcdir}/${slcdate}.rslc.mli.par | awk '{print $2}'` 
    reducfac=`echo $wid | awk '{if(int($1/1000) > 1) print int($1/1000); else print 1}'` 
    raspwr ${rslcdir}/${slcdate}.rslc.mli $wid - - $reducfac $reducfac - - - ${rslcdir}/${slcdate}.rslc.mli.ras 0 > /dev/null 2>&1 #>> $logfile 
  fi
}; export -f multilookRSLC

########################################################################
########################################################################
# remove all original unzipped slc files from the SLC/yyyymmdd directory
function rmTmpSLC(){
  local slcdir=$1
  # Remove the original files associated the SAFE files 
  rm -f $slcdir/????????T*.IW*.slc    $slcdir/????????T*.IW*.slc.par    $slcdir/????????T*.IW*.slc.TOPS_par
  # Remove the intermediate files used to concatenate more than 2 slice products  
  rm -f $slcdir/????????.[0-30].*.slc $slcdir/????????.[0-30].*.slc.par $slcdir/????????.[0-30].*.slc.TOPS_par
}; export -f rmTmpSLC

function check_master(){
  #######################################################
  # check_master list_of_zipfiles master_date
  #######################################################
  local zipfile=$1
  local master_date=$2
  if [ -z "$master_date" ] ; then
    # We choose the first image as master, because master_date is not available
    slclist=( `cat ${zipfile}` ) ; sortarray_slc_yyyymmdd slclist[@] SLC/ ; 
    master_date=`basename $(echo ${slclist[0]} ) | awk '{print substr($1,18,8)}'`
    echo "$master_date" > master_date.txt
  else
    echo " Master date introduced by user: $master_date "
  fi
}; export -f check_master # End function definition


function check_EmptySubswath(){
  local image1=$1  
  if [ -z "$2" ]; then local IW1=1; else local IW1=$2; fi 
  if [ -z "$3" ]; then local IW2=3; else local IW2=$3; fi 
  rm -f temp
  for i in `seq ${IW1} ${IW2}`; do
    MinBurstIW=`awk '{print $1}' ${image1}_IW${i}.commonburst | sort -n | head -1`;
    MaxBurstIW=`awk '{print $2}' ${image1}_IW${i}.commonburst | sort -n | tail -1`; 
    echo $MinBurstIW $MaxBurstIW | awk '{if($1!=0 || $2!=0) print "'${i}'"}' >> temp
  done
  echo `sort -n temp | head -1` `sort -n temp | tail -1`
}; export -f check_EmptySubswath # End function definition

########################################################################
########################################################################
# Move files from a cropped version to normal version (Merging step)
function mk_masterRSLClinks(){
  local iput=$1
  local oput=$2
  if [ -z "$3" ]; then local IW1=1; else local IW1=$3; fi 
  if [ -z "$4" ]; then local IW2=3; else local IW2=$4; fi   
  ln -s ${iput}.slc     ${oput}.rslc 
  ln -s ${iput}.slc.par ${oput}.rslc.par 
  for i in `seq ${IW1} ${IW2}`; do
    ln -s ${iput}.IW${i}.slc           ${oput}.IW${i}.rslc 
    ln -s ${iput}.IW${i}.slc.par       ${oput}.IW${i}.rslc.par 
    ln -s ${iput}.IW${i}.slc.TOPS_par  ${oput}.IW${i}.rslc.TOPS_par
  done
}; export -f mk_masterRSLClinks # End function definition

function mk_RSLClinks(){
  local iput=$1
  local iend=$2
  local oput=$3
  local oend=$4
  if [ -z "$5" ]; then local IW1=1; else local IW1=$5; fi 
  if [ -z "$6" ]; then local IW2=3; else local IW2=$6; fi   
  ln -s ${iput}.${iend}     ${oput}.${oend} 
  ln -s ${iput}.${iend}.par ${oput}.${oend}.par 
  for i in `seq ${IW1} ${IW2}`; do
    ln -s ${iput}.IW${i}.${iend}           ${oput}.IW${i}.${oend} 
    ln -s ${iput}.IW${i}.${iend}.par       ${oput}.IW${i}.${oend}.par 
    ln -s ${iput}.IW${i}.${iend}.TOPS_par  ${oput}.IW${i}.${oend}.TOPS_par
  done
}; export -f mk_RSLClinks # End function definition

function mk_RMLISLClinks(){
  #mk_RMLISLClinks ../SLC/${master} slc.mli ${master} rmli ${IW1} ${IW3}
  local iput=$1
  local iend=$2
  local oput=$3
  local oend=$4
  if [ -z "$5" ]; then local IW1=1; else local IW1=$5; fi 
  if [ -z "$6" ]; then local IW2=3; else local IW2=$6; fi   
  ln -s ${iput}.${iend}     ${oput}.${oend} 
  ln -s ${iput}.${iend}.par ${oput}.${oend}.par 
}; export -f mk_RMLISLClinks # End function definition

function mk_updatelt(){
  # Update the lookup table!
  # This section is computed to update: lookup table refinement
  # determine range and azimuth corrections for lookup table (in mli pixels)
  
  local inlt=$1
  local outlt=$2
  local filein=$3 # ${master2slave}.fine.off
  local master=$4 
  local master2slave=$5
  local width=$6
  local rlks=$7
  local azlks=$8
  
  dr=`awk '$1 == "range_offset_polynomial:" {print $2}' ${filein}`
  dr_mli=`echo "$dr" "$rlks" | awk '{printf "%f", $1/$2}'`
  daz=`awk '$1 == "azimuth_offset_polynomial:" {print $2}' ${filein}`
  daz_mli=`echo "$daz" "$azlks" | awk '{printf "%f", $1/$2}'`
  create_diff_par SLC/${master}.slc.mli.par SLC/${master}.slc.mli.par ${master2slave}.diff_par 1 0 &> /dev/null
  set_value ${master2slave}.diff_par ${master2slave}.diff_par "range_offset_polynomial"   "$dr_mli   0.0000e+00   0.0000e+00   0.0000e+00   0.0000e+00   0.0000e+00" &> /dev/null
  set_value ${master2slave}.diff_par ${master2slave}.diff_par "azimuth_offset_polynomial" "$daz_mli   0.0000e+00   0.0000e+00   0.0000e+00   0.0000e+00   0.0000e+00" &> /dev/null
  cp ${master2slave}.diff_par ${master2slave}.diff_par.1
  mv ${inlt} RSLC/${master}.slc.mli.lt.tmp.1
  gc_map_fine RSLC/${master}.slc.mli.lt.tmp.1 $width ${master2slave}.diff_par ${outlt} 1 &> /dev/null
  mv RSLC/${master}.slc.mli.lt.tmp.1 ${inlt} 
}; export -f mk_updatelt # End function definition

function plotifg(){
  local master=$1
  local slave=$2
  local inputgrd=geocode/${master}_${slave}.geo.100m.cm.grd
  gmt gmtset HEADER_FONT_SIZE 10 HEADER_OFFSET 0.15 ANNOT_FONT_SIZE_SECONDARY 10 PLOT_DEGREE_FORMAT DF BASEMAP_TYPE plain CHAR_ENCODING Standard+ \
       COLOR_MODEL RGB DOTS_PR_INCH 1200 INPUT_DATE_FORMAT yyyy/mm/dd PLOT_DATE_FORMAT o \
       TIME_FORMAT_PRIMARY abbreviated ANNOT_FONT_SIZE_PRIMARY 7p LABEL_FONT_SIZE 7p \
       TICK_LENGTH 0.1c ANNOT_FONT_SIZE_PRIMARY 7p ANNOT_OFFSET_PRIMARY 0.1c \
       LABEL_FONT_SIZE 7p LABEL_OFFSET -0.05c  
  grd2cpt $inputgrd -D -I -E100 -C${LiCSARpath}/misc/mypolar.cpt > color.cpt
  grdimage ${inputgrd} -JM16 -R${inputgrd} -Ccolor.cpt -Q -P -K > geocode/${master}_${slave}.geo.cm.ps
  psscale -D8/1.0/6/0.3h -Ccolor.cpt -B4:"Relative motion towards/away satellite  [cm]": -O >> geocode/${master}_${slave}.geo.cm.ps
  ps2raster geocode/${master}_${slave}.geo.cm.ps -E600 -Tg -W+k+t"S1 ${master}_${slave}"+l256/-1
  convert geocode/${master}_${slave}.geo.cm.png -transparent white geocode/${master}_${slave}.geo.cm.png
  rm -f color.cpt geocode/${master}_${slave}.geo.cm.ps 
}; export -f plotifg

function plotdem(){
  local master=$1
  local slave=$2
  local inputgrd=geocode/${master}_${slave}.geo.100m.dem.grd
  gmt gmtset HEADER_FONT_SIZE 10 HEADER_OFFSET 0.15 ANNOT_FONT_SIZE_SECONDARY 10 PLOT_DEGREE_FORMAT DF BASEMAP_TYPE plain CHAR_ENCODING Standard+ \
       COLOR_MODEL RGB DOTS_PR_INCH 1200 INPUT_DATE_FORMAT yyyy/mm/dd PLOT_DATE_FORMAT o \
       TIME_FORMAT_PRIMARY abbreviated ANNOT_FONT_SIZE_PRIMARY 7p LABEL_FONT_SIZE 7p \
       TICK_LENGTH 0.1c ANNOT_FONT_SIZE_PRIMARY 7p ANNOT_OFFSET_PRIMARY 0.1c \
       LABEL_FONT_SIZE 7p LABEL_OFFSET -0.05c  
  grdimage ${inputgrd} -JM16 -R${inputgrd} -C${LiCSARpath}/misc/nicetopo.cpt -Q -P -K > geocode/${master}_${slave}.geo.dem.ps
  psscale -D8/1.0/6/0.3h -C/nfs/see-fs-01_users/earpjg/pablo/CommonData/CPTs/nicetopo.cpt -B1000:"Topography [m]": -O >> geocode/${master}_${slave}.geo.dem.ps
  ps2raster geocode/${master}_${slave}.geo.dem.ps -E600 -Tg -W+k+t"S1 ${master}_${slave}"+l256/-1
  convert geocode/${master}_${slave}.geo.dem.png -transparent white geocode/${master}_${slave}.geo.dem.png
  rm -f color.cpt geocode/${master}_${slave}.geo.dem.ps
}; export -f plotdem

function LiCSAR_plotdem(){
  local master=$1
  local inputgrd=geocode/${master}.geo.dem.grd
  gmt gmtset HEADER_FONT_SIZE 10 HEADER_OFFSET 0.15 ANNOT_FONT_SIZE_SECONDARY 10 PLOT_DEGREE_FORMAT DF BASEMAP_TYPE plain CHAR_ENCODING Standard+ \
       COLOR_MODEL RGB DOTS_PR_INCH 1200 INPUT_DATE_FORMAT yyyy/mm/dd PLOT_DATE_FORMAT o \
       TIME_FORMAT_PRIMARY abbreviated ANNOT_FONT_SIZE_PRIMARY 7p LABEL_FONT_SIZE 7p \
       TICK_LENGTH 0.1c ANNOT_FONT_SIZE_PRIMARY 7p ANNOT_OFFSET_PRIMARY 0.1c \
       LABEL_FONT_SIZE 7p LABEL_OFFSET -0.05c  
  grdimage ${inputgrd} -JM16 -R${inputgrd} -C${LiCSARpath}/misc/nicetopo.cpt -Q -P -K > geocode/${master}.geo.dem.ps
  psscale -D8/1.0/6/0.3h -C/nfs/see-fs-01_users/earpjg/pablo/CommonData/CPTs/nicetopo.cpt -B1000:"Topography [m]": -O >> geocode/${master}.geo.dem.ps
  ps2raster geocode/${master}.geo.dem.ps -E600 -Tg -W+k+t"S1 ${master}"+l256/-1
  convert geocode/${master}.geo.dem.png -transparent white geocode/${master}.geo.dem.png
  rm -f color.cpt geocode/${master}.geo.dem.ps
}; export -f LiCSAR_plotdem


function LiCSAR_plotifg(){
  local master=$1
  local slave=$2
  local workdir=$3
  local extension=$4
  local transparencyon=$5
  local inputgrd=${workdir}/${master}_${slave}.${extension}.grd
  rm -f .gmtdefaults4 .gmtcommands4
  gmt gmtset HEADER_FONT_SIZE 10 HEADER_OFFSET 0.15 ANNOT_FONT_SIZE_SECONDARY 10 PLOT_DEGREE_FORMAT DF BASEMAP_TYPE plain CHAR_ENCODING Standard+ \
         COLOR_MODEL RGB DOTS_PR_INCH 1200 INPUT_DATE_FORMAT yyyy/mm/dd PLOT_DATE_FORMAT o \
         TIME_FORMAT_PRIMARY abbreviated ANNOT_FONT_SIZE_PRIMARY 7p LABEL_FONT_SIZE 7p \
         TICK_LENGTH 0.1c ANNOT_FONT_SIZE_PRIMARY 7p ANNOT_OFFSET_PRIMARY 0.1c \
         LABEL_FONT_SIZE 7p LABEL_OFFSET -0.05c # COLOR_BACKGROUND 50/50/50 COLOR_NAN 50/50/50 
  grd2cpt $inputgrd -D -I -E100 -C${LiCSARpath}/misc/mypolar.cpt > color.cpt
  maxvalue=`grdinfo -C ${inputgrd} | awk '{print $7}'`
  minvalue=`grdinfo -C ${inputgrd} | awk '{print $6}'`
  stepticks=`echo $maxvalue $minvalue | awk '{print int(($1-$2)/8)}'`
  grdmath ${inputgrd} 0.0 NAN = temp.grd
  grdimage temp.grd -JM16 -R${inputgrd} -Ccolor.cpt -Q -P -K > ${workdir}/${master}_${slave}.${extension}.ps
  psscale -D8/1.0/6/0.3h -Ccolor.cpt -B${stepticks}:"Relative motion towards/away satellite  [cm]": -O >> ${workdir}/${master}_${slave}.${extension}.ps
  ps2raster ${workdir}/${master}_${slave}.${extension}.ps -E600 -Tg -W+k+t"S1 ${master}_${slave}"+l256/-1
  if [ "$transparencyon" == "1" ]; then
    #convert ${workdir}/${master}_${slave}.${extension}.png -fuzz 1% -transparent "rgb(50,50,50)" ${workdir}/${master}_${slave}.${extension}.png
    convert ${workdir}/${master}_${slave}.${extension}.png -transparent white ${workdir}/${master}_${slave}.${extension}.png
  fi
  rm -f color.cpt geocode/${master}_${slave}.${extension}.ps temp.grd
}; export -f LiCSAR_plotifg

function mk_envihdr(){
  local width_dem=$1 
  local length_dem=$2
  local lon=$3
  local lat=$4
  local lonstep=$5
  local latstep=$6
  
  echo "ENVI description = { Registration Result. Method1st degree Polynomial w/ nearest neighbor [Wed Dec 20 23:59:19 1995] }"
  echo "samples = $width_dem"
  echo "lines   = $length_dem"
  echo "bands   = 1"
  echo "header offset = 0"
  echo "file type = ENVI Standard"
  echo "data type = 4"
  echo "interleave = bsq"
  echo "byte order = 1"
  echo "map info = {Geographic Lat/Lon, 1, 1, ${lon}, ${lat}, ${lonstep}, ${latstep},  WGS-84}"
  echo "coordinate system string = {GEOGCS["GCS_WGS_1984",DATUM["D_WGS_1984",SPHEROID["WGS_1984",6378137,298.257223563]],PRIMEM["Greenwich",0],UNIT["Degree",0.017453292519943295]]}"
  echo "band names = {Band 1}"
}; export -f mk_envihdr

function mk_squint_TOPSmosaic(){
  local master=$1 # Mosaic 
  local outfile=$2
  local rlks=$3
  local azlks=$4
  local wid=$5  
  if [ -z "$6" ]; then local IW1=1; else local IW1=$6; fi 
  if [ -z "$7" ]; then local IW2=3; else local IW2=$7; fi  
  
  # Loop over all subswaths of the master
  for i in `seq ${IW1} ${IW2}`; do
    # Simulate squint angles and write a complex file 
    # Read PRF and steering rate
    #PRF=`awk '$1 == "prf:" {print $2}' SLC/${master}.IW${i}.slc.par` # in Hz
    #K=`awk '$1 == "az_steering_rate:" {print $2}' SLC/${master}.IW${i}.slc.TOPS_par` # in deg/sec
    
    length_file=`awk '$1 == "azimuth_lines:" {print $2}' SLC/${master}.IW${i}.slc.par` # in pixels
    azimuth_pixel_spacing=`awk '$1 == "azimuth_pixel_spacing:" {print $2}' SLC/${master}.IW${i}.slc.par` # in meters
    center_range_distance=`awk '$1 == "center_range_slc:" {print $2}' SLC/${master}.IW${i}.slc.par` # in meters
    nlines_per_burst=`awk '$1 == "lines_per_burst:" {print $2}' SLC/${master}.IW${i}.slc.TOPS_par` # in pixels
    width_burst=`awk '$1 == "range_samples:" {print $2}' SLC/${master}.IW${i}.slc.par` # in pixels
    nbursts=`awk '$1 == "number_of_bursts:" {print $2}' SLC/${master}.IW${i}.slc.TOPS_par` # in deg/sec
    
    echo "Creating Squint angles for master geometry: SLC/${master}.IW${i}"
    echo "  mk_beta_mat_alternativemethod length_file azimuth_pixel_spacing center_range_slc nlines_per_burst width_burst nbursts output "
    echo "  mk_beta_mat_alternativemethod $length_file $azimuth_pixel_spacing $center_range_distance $nlines_per_burst $width_burst $nbursts SLC/${master}.IW${i}.slc.beta "
    # Call function that calculates a large matrix in complex format with phase equals to the squint angle
    #mk_beta_mat SLC/${master}.IW${i}.slc $length_file $PRF $K $nlines_per_burst $width_burst SLC/${master}.IW${i}.slc.beta 
    mk_beta_mat_alternativemethod $length_file $azimuth_pixel_spacing $center_range_distance $nlines_per_burst $width_burst $nbursts SLC/${master}.IW${i}.slc.beta 

  done
  
  # Create fake SLC_tab files
  createBetaSLCtab SLC/${master} slc ${IW1} ${IW2} > BetaSLC_tab
  
  # Merge all simulated complex files with the squint angles
  echo "Mosaicing bursted squint angles for master geometry: SLC/${master}"
  SLC_mosaic_S1_TOPS BetaSLC_tab ${outfile} ${outfile}.par $rlks $azlks 
  
  # Multilook the results
  echo "Multilook mosaic squint angle matrix for master geometry: SLC/${master}"
  create_offset SLC/${master}.slc.par SLC/${master}.slc.par offsetfile.off 1 1 1 0  
  multi_cpx ${outfile} offsetfile.off ${outfile}.ml offsetfile.off $rlks $azlks

  # Extract the phase from the multilooked complex files
  echo "Extract phase values from complex squint angle matrix file for master geometry: SLC/${master}"
  cpx_to_real ${outfile}.ml ${outfile} $wid 4

  # Cleaning
  for i in `seq ${IW1} ${IW2}`; do rm -f SLC/${master}.IW${i}.slc.beta ; done
  rm -f BetaSLC_tab ${outfile} ${outfile}.par offsetfile.off ${outfile}.ml 

}; export -f mk_squint_TOPSmosaic

########################################################################
########################################################################
# 
function createBetaSLCtab(){
  local inimage=$1
  local ending=$2
  if [ -z "$3" ]; then local IW1=1; else local IW1=$3; fi 
  if [ -z "$4" ]; then local IW2=3; else local IW2=$4; fi   
  for i in `seq ${IW1} ${IW2}`; do
    echo "${inimage}.IW${i}.${ending}.beta ${inimage}.IW${i}.${ending}.par ${inimage}.IW${i}.${ending}.TOPS_par"
  done
}; export -f createBetaSLCtab

function mk_beta_mat_alternativemethod(){
  # mk_beta_mat_alternativemethod $length_file $azimuth_pixel_spacing $center_range_slc $nlines_per_burst $width_burst $nbursts SLC/${master}.IW${i}.slc.beta 
  local length=$1
  local azimuth_pixel_spacing=$2
  local center_range_distance=$3
  local nlines_per_burst=$4
  local width_burst=$5
  local nbursts=$6
  local outfile=$7
  
cat << __EOF__ >> matlab_func.m  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Script to compute squint angle matrix (output in radians)
%
% Author: Pablo J. Gonzalez
% Date: 2015/12/22
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
SLCbeta = ones(${length},${width_burst}); 
SLCbeta = complex(SLCbeta,0);
%beta_range = rad2deg ( ( (${azimuth_pixel_spacing}*${nlines_per_burst}) ) / ${center_range_distance} );
beta_range = (${azimuth_pixel_spacing}*${nlines_per_burst}) / ${center_range_distance} ;
betavec = -(beta_range/2):beta_range/(${nlines_per_burst}-1):(beta_range/2);
betamat = repmat(betavec,${width_burst},1)';
betaMAT = repmat(betamat,${nbursts},1);
SLCbeta=SLCbeta.*exp(1i*betaMAT);
fwritebk(SLCbeta,'${outfile}','cpxfloat32','b');
__EOF__

matlab -nojvm -nosplash < matlab_func.m &> /dev/null 

rm -f matlab_func.m 
}; export -f mk_beta_mat_alternativemethod


function mk_geom_coreg_TOPS(){
  local master=$1
  local slave=$2
  local demdir=$3
  if [ -z "$4" ]; then local IW1=1; else local IW1=$4; fi 
  if [ -z "$5" ]; then local IW2=3; else local IW2=$5; fi 
  if [ -z "$6" ]; then local logfile=dummy.log; else local logfile=$6; fi 
  
  createSLCtab SLC/${master}/${master} slc  ${IW1} ${IW3} > SLC1_tab
  createSLCtab SLC/${slave}/${slave}  slc  ${IW1} ${IW3} > SLC2_tab
  createSLCtab RSLC/${slave}/${slave} rslc ${IW1} ${IW3} > RSLC2_tab
  
  # Coregister images using geometric corregistration
  rdc_trans SLC/${master}/${master}.slc.mli.par ${demdir}/${master}.hgt SLC/${slave}/${slave}.slc.mli.par RSLC/${slave}/${master}_${slave}.slc.mli.lt >> $logfile
  SLC_interp_lt_S1_TOPS SLC2_tab SLC/${slave}/${slave}.slc.par SLC1_tab SLC/${master}/${master}.slc.par RSLC/${slave}/${master}_${slave}.slc.mli.lt SLC/${master}/${master}.slc.mli.par SLC/${slave}/${slave}.slc.mli.par - RSLC2_tab RSLC/${slave}/${slave}.rslc RSLC/${slave}/${slave}.rslc.par >> $logfile

}; export -f mk_geom_coreg_TOPS


function mk_geom_coreg_TOPS_recropmaster(){
  local master=$1
  local slave=$2
  local demdir=$3
  if [ -z "$4" ]; then local IW1=1; else local IW1=$4; fi 
  if [ -z "$5" ]; then local IW2=3; else local IW2=$5; fi 
  if [ -z "$6" ]; then local logfile=dummy.log; else local logfile=$6; fi 
  
  createSLCtab RSLC/${slave}/${master} slc  ${IW1} ${IW3} > SLC1_tab
  createSLCtab SLC/${slave}/${slave}   slc  ${IW1} ${IW3} > SLC2_tab
  createSLCtab RSLC/${slave}/${slave}  rslc ${IW1} ${IW3} > RSLC2_tab
  
  # Coregister images using geometric corregistration
  rdc_trans RSLC/${slave}/${master}.slc.mli.par ${demdir}/${master}.hgt SLC/${slave}/${slave}.slc.mli.par RSLC/${slave}/${master}_${slave}.slc.mli.lt >> $logfile
  SLC_interp_lt_S1_TOPS SLC2_tab SLC/${slave}/${slave}.slc.par SLC1_tab RSLC/${slave}/${master}.slc.par RSLC/${slave}/${master}_${slave}.slc.mli.lt RSLC/${slave}/${master}.slc.mli.par SLC/${slave}/${slave}.slc.mli.par - RSLC2_tab RSLC/${slave}/${slave}.rslc RSLC/${slave}/${slave}.rslc.par >> $logfile

}; export -f mk_geom_coreg_TOPS_recropmaster


function mk_coreg_TOPS_specdiv(){
  local master=$1
  local slave=$2
  local width=$3
  local rlks=$4
  local azlks=$5
  if [ -z "$6" ]; then local logfile=dummy.log; else local logfile=$6; fi 
  if [ -z "$7" ]; then local IW1=1; else local IW1=$7; fi 
  if [ -z "$8" ]; then local IW3=3; else local IW3=$8; fi
    
  createSLCtab SLC/${master}/${master}  slc  ${IW1} ${IW3} > SLC1_tab
  createSLCtab SLC/${slave}/${slave}    slc  ${IW1} ${IW3} > SLC2_tab
  createSLCtab RSLC/${slave}/${slave}   rslc ${IW1} ${IW3} > RSLC2_tab
  
  # Create the offset file
  create_offset SLC/${master}/${master}.slc.par SLC/${slave}/${slave}.slc.par RSLC/${slave}/${master}_${slave}.off 1 ${rlks} ${azlks} 0 >> $logfile
  
  # If a closer in time slave exists, use it to compute spectral diversity
  if [ -z "$9" ]; then 
    echo "   Resampling using spectral diversity master-slave [${master} - ${slave}]"
  else
    echo "   Resampling using spectral diversity slave2-slave [${slave3} - ${slave}]"  
    local slave3=$9
    createSLCtab RSLC/${slave3}/${slave3} rslc ${IW1} ${IW3} > RSLC3_tab
  fi  
  
  # Refine twice offset with spectral diversity
  if [ ! -f RSLC3_tab ]; then 
    S1_coreg_overlap SLC1_tab RSLC2_tab RSLC/${slave}/${master}_${slave} RSLC/${slave}/${master}_${slave}.off RSLC/${slave}/${master}_${slave}.off.refine1 0.8 0.01 0.8 1 >> $logfile #2>&1 $logfile  
    SLC_interp_lt_S1_TOPS SLC2_tab SLC/${slave}/${slave}.slc.par SLC1_tab SLC/${master}/${master}.slc.par RSLC/${slave}/${master}_${slave}.slc.mli.lt SLC/${master}/${master}.slc.mli.par SLC/${slave}/${slave}.slc.mli.par RSLC/${slave}/${master}_${slave}.off.refine1 RSLC2_tab RSLC/${slave}/${slave}.rslc RSLC/${slave}/${slave}.rslc.par >> $logfile
    S1_coreg_overlap SLC1_tab RSLC2_tab RSLC/${slave}/${master}_${slave} RSLC/${slave}/${master}_${slave}.off.refine1 RSLC/${slave}/${master}_${slave}.off.refine2 0.8 0.01 0.8 1 >> $logfile #2>&1 $logfile
    SLC_interp_lt_S1_TOPS SLC2_tab SLC/${slave}/${slave}.slc.par SLC1_tab SLC/${master}/${master}.slc.par RSLC/${slave}/${master}_${slave}.slc.mli.lt SLC/${master}/${master}.slc.mli.par SLC/${slave}/${slave}.slc.mli.par RSLC/${slave}/${master}_${slave}.off.refine2 RSLC2_tab RSLC/${slave}/${slave}.rslc RSLC/${slave}/${slave}.rslc.par >> $logfile    
    rm -f RSLC/${slave}/${master}_${slave}.IW*.*.diff* RSLC/${slave}/${master}_${slave}.IW*.*.int* RSLC/${slave}/${master}_${slave}.IW*.*.off* 
  else
    S1_coreg_overlap SLC1_tab RSLC2_tab RSLC/${slave}/${master}_${slave} RSLC/${slave}/${master}_${slave}.off RSLC/${slave}/${master}_${slave}.off.refine1 0.8 0.01 0.8 1 RSLC3_tab >> $logfile #2>&1 $logfile
    SLC_interp_lt_S1_TOPS SLC2_tab SLC/${slave}/${slave}.slc.par SLC1_tab SLC/${master}/${master}.slc.par RSLC/${slave}/${master}_${slave}.slc.mli.lt SLC/${master}/${master}.slc.mli.par SLC/${slave}/${slave}.slc.mli.par RSLC/${slave}/${master}_${slave}.off.refine1 RSLC2_tab RSLC/${slave}/${slave}.rslc RSLC/${slave}/${slave}.rslc.par >> $logfile
    S1_coreg_overlap SLC1_tab RSLC2_tab RSLC/${slave}/${master}_${slave} RSLC/${slave}/${master}_${slave}.off.refine1 RSLC/${slave}/${master}_${slave}.off.refine2 0.8 0.01 0.8 1 RSLC3_tab >> $logfile #2>&1 $logfile
    SLC_interp_lt_S1_TOPS SLC2_tab SLC/${slave}/${slave}.slc.par SLC1_tab SLC/${master}/${master}.slc.par RSLC/${slave}/${master}_${slave}.slc.mli.lt SLC/${master}/${master}.slc.mli.par SLC/${slave}/${slave}.slc.mli.par RSLC/${slave}/${master}_${slave}.off.refine2 RSLC2_tab RSLC/${slave}/${slave}.rslc RSLC/${slave}/${slave}.rslc.par >> $logfile
    rm -f RSLC/${slave}/${master}_${slave}.IW*.*.diff* RSLC/${slave}/${master}_${slave}.IW*.*.int* RSLC/${slave}/${master}_${slave}.IW*.*.off* 
  fi
}; export -f mk_coreg_TOPS_specdiv

function mk_coreg_TOPS_specdiv_crop(){
  local master=$1
  local slave=$2
  local width=$3
  local rlks=$4
  local azlks=$5
  if [ -z "$6" ]; then local logfile=dummy.log; else local logfile=$6; fi 
  if [ -z "$7" ]; then local IW1=1; else local IW1=$7; fi 
  if [ -z "$8" ]; then local IW3=3; else local IW3=$8; fi
    
  createSLCtab RSLC/${slave}/${master}  slc  ${IW1} ${IW3} > SLC1_tab
  createSLCtab SLC/${slave}/${slave}    slc  ${IW1} ${IW3} > SLC2_tab
  createSLCtab RSLC/${slave}/${slave}   rslc ${IW1} ${IW3} > RSLC2_tab
  
  # Create the offset file
  create_offset RSLC/${slave}/${master}.slc.par SLC/${slave}/${slave}.slc.par RSLC/${slave}/${master}_${slave}.off 1 ${rlks} ${azlks} 0 >> $logfile
  
#   # If a closer in time slave exists, use it to compute spectral diversity
#   if [ -z "$9" ]; then 
#     echo "   Resampling using spectral diversity master-slave [${master} - ${slave}]"
#   else
#     echo "   Resampling using spectral diversity slave2-slave [${slave3} - ${slave}]"  
#     local slave3=$9
#     createSLCtab RSLC/${slave3}/${slave3} rslc ${IW1} ${IW3} > RSLC3_tab
#   fi  
  
  # Refine twice offset with spectral diversity
#   if [ ! -f RSLC3_tab ]; then 
    S1_coreg_overlap SLC1_tab RSLC2_tab RSLC/${slave}/${master}_${slave} RSLC/${slave}/${master}_${slave}.off RSLC/${slave}/${master}_${slave}.off.refine1 0.8 0.01 0.8 1 >> $logfile #2>&1 $logfile  
    SLC_interp_lt_S1_TOPS SLC2_tab SLC/${slave}/${slave}.slc.par SLC1_tab SLC/${master}/${master}.slc.par RSLC/${slave}/${master}_${slave}.slc.mli.lt RSLC/${slave}/${master}.slc.mli.par SLC/${slave}/${slave}.slc.mli.par RSLC/${slave}/${master}_${slave}.off.refine1 RSLC2_tab RSLC/${slave}/${slave}.rslc RSLC/${slave}/${slave}.rslc.par >> $logfile
    S1_coreg_overlap SLC1_tab RSLC2_tab RSLC/${slave}/${master}_${slave} RSLC/${slave}/${master}_${slave}.off.refine1 RSLC/${slave}/${master}_${slave}.off.refine2 0.8 0.01 0.8 1 >> $logfile #2>&1 $logfile
    SLC_interp_lt_S1_TOPS SLC2_tab SLC/${slave}/${slave}.slc.par SLC1_tab SLC/${master}/${master}.slc.par RSLC/${slave}/${master}_${slave}.slc.mli.lt RSLC/${slave}/${master}.slc.mli.par SLC/${slave}/${slave}.slc.mli.par RSLC/${slave}/${master}_${slave}.off.refine2 RSLC2_tab RSLC/${slave}/${slave}.rslc RSLC/${slave}/${slave}.rslc.par >> $logfile    
    rm -f RSLC/${slave}/${master}_${slave}.IW*.*.diff* RSLC/${slave}/${master}_${slave}.IW*.*.int* RSLC/${slave}/${master}_${slave}.IW*.*.off* 
#   else
#     S1_coreg_overlap SLC1_tab RSLC2_tab RSLC/${slave}/${master}_${slave} RSLC/${slave}/${master}_${slave}.off RSLC/${slave}/${master}_${slave}.off.refine1 0.8 0.01 0.8 1 RSLC3_tab >> $logfile #2>&1 $logfile
#     SLC_interp_lt_S1_TOPS SLC2_tab SLC/${slave}/${slave}.slc.par SLC1_tab SLC/${master}/${master}.slc.par RSLC/${slave}/${master}_${slave}.slc.mli.lt SLC/${master}/${master}.slc.mli.par SLC/${slave}/${slave}.slc.mli.par RSLC/${slave}/${master}_${slave}.off.refine1 RSLC2_tab RSLC/${slave}/${slave}.rslc RSLC/${slave}/${slave}.rslc.par >> $logfile
#     S1_coreg_overlap SLC1_tab RSLC2_tab RSLC/${slave}/${master}_${slave} RSLC/${slave}/${master}_${slave}.off.refine1 RSLC/${slave}/${master}_${slave}.off.refine2 0.8 0.01 0.8 1 RSLC3_tab >> $logfile #2>&1 $logfile
#     SLC_interp_lt_S1_TOPS SLC2_tab SLC/${slave}/${slave}.slc.par SLC1_tab SLC/${master}/${master}.slc.par RSLC/${slave}/${master}_${slave}.slc.mli.lt SLC/${master}/${master}.slc.mli.par SLC/${slave}/${slave}.slc.mli.par RSLC/${slave}/${master}_${slave}.off.refine2 RSLC2_tab RSLC/${slave}/${slave}.rslc RSLC/${slave}/${slave}.rslc.par >> $logfile
#     rm -f RSLC/${slave}/${master}_${slave}.IW*.*.diff* RSLC/${slave}/${master}_${slave}.IW*.*.int* RSLC/${slave}/${master}_${slave}.IW*.*.off* 
#   fi
}; export -f mk_coreg_TOPS_specdiv_crop


function deramp_ifg {
 path=$1
 if [ -z $1 ]; then echo "parameter is full path to the ifg geotiff in LiCSAR_public"; 
  else
    python3 -c "import LiCSAR_lib.unwrp_multiscale as unw; unw.deramp_ifg_tif('"$path"')"
  fi
}; export -f deramp_ifg

function datediff {
# very useful function to find the difference (in days) between two dates
A=$1
B=$2
#B=`date +%Y-%m-%d`
#A=${1:0:4}-${1:4:2}-${1:6:2}
#B=${2:0:4}-${2:4:2}-${2:6:2}
sec=86400
dte1=$(date --utc --date "$A" +%s)
dte2=$(date --utc --date "$B" +%s)
diffSec=$((dte2-dte1))
if ((diffSec < 0)); then abs=-1; else abs=1; fi
ROZDIL=$((diffSec/sec*abs))
echo $ROZDIL
}; export -f datediff
