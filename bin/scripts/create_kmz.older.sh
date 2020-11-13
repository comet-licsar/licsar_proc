#!/bin/bash
#this will force regenerate full previews
overwrite=1
update_previews=0

if [ $overwrite -eq 1 ]; then
 update_previews=1
fi

if [ -z $1 ]; then echo "parameter is full path to the folder containing geotiffs";
echo "e.g. /gws/nopw/j04/nceo_geohazards_vol1/public/LiCSAR_products/58/058A_05279_131311/interferograms/20190317_20190323"
echo "if second (optional) parameter is a frame name, it will use it in the output filename (e.g. 058A_05279_131311_20200505_20200517.kmz)"
exit;
fi

if [ ! -z $2 ]; then opt=$2; else opt=0; fi

pair=`basename $1`
if [ ! `echo $pair | wc -m` -eq 18 ]; then 
 echo "wrong path - it should end on pairname, e.g. 20190101_20190107"
 exit
fi

if [ ! -f $1/$pair.geo.unw.png ] || [ ! -f $1/$pair.geo.diff.png ] || [ ! -f $1/$pair.geo.unw.tif ]; then
 echo "some of the geocoded product does not exist here"
 exit
fi

### let's make it
echo "generating kmz file for pair "$pair
cd $1
unw_tif=$pair.geo.unw.tif
unw_bmp=$pair.geo.unw.png
#ifg_tif=$pair.geo.diff_pha.tif
#ifg_bmp=$pair.geo.diff.png
#ifg_unfiltered_bmp=$pair.geo.diff_unfiltered.png
coh_bmp=$pair.geo.cc.png
ifg_bmp=$pair.geo.diff.full.png
ifg_unfiltered_bmp=$pair.geo.diff_unfiltered.full.png

if [ -f $pair.geo.cc.full.png ]; then
 coh_bmp=$pair.geo.cc.full.png
fi
if [ -f $pair.geo.unw.full.png ]; then
 unw_bmp=$pair.geo.unw.full.png;
fi


do_ifg=0
do_ifgu=0
do_unw=0
if [ $overwrite -eq 1 ]; then
 do_ifg=1
 do_ifgu=1
 do_unw=1
 echo "generating fullscale previews"
 #doing also unw for full scale then
 unw_bmp=$pair.geo.unw.full.png
 if [ ! -f $pair.geo.cc.full.png ]; then
   #regenerate full coh image..
   gmt makecpt -Cgray -T0/255/1 >cc.cpt
   gmt grdimage $pair.geo.cc.tif -Ccc.cpt -JM1 -nn+t0.1 -Q -Abbb.png
   coh_bmp=$pair.geo.cc.full.png
   convert -transparent black bbb.png PNG8:$coh_bmp
   rm bbb.png cc.cpt
 fi
else
 if [ ! -f $ifg_bmp ]; then
  do_ifg=1
 fi
 if [ ! -f $ifg_unfiltered_bmp ]; then
  do_ifgu=1
 fi
 if [ ! -f $unw_bmp ]; then
   do_unw=1
 fi
fi

if [ $do_ifg -eq 1 ]; then
 echo "picturing filtered ifg"
 #first doing filtered ifg
 gmt grdimage $pair.geo.diff_pha.tif -C$LiCSARpath/misc/pha.cpt -JM1 -nn+t0.1 -Q -A$ifg_bmp
 #to flatten it (and fix transparency...sometimes needed..):
 convert $ifg_bmp -transparent black PNG8:temptemp.png
 mv temptemp.png $ifg_bmp
fi

if [ $do_ifgu -eq 1 ]; then
 #perhaps we are to add unfiltered one:
  if [ -f $pair.geo.diff_unfiltered_pha.tif ]; then
    echo "picturing unfiltered ifg"
    gmt grdimage $pair.geo.diff_unfiltered_pha.tif -C$LiCSARpath/misc/pha.cpt -Q -JM1 -nn+t0.1 -A$ifg_unfiltered_bmp
    convert $ifg_unfiltered_bmp -transparent black PNG8:temptemp.png
    mv temptemp.png $ifg_unfiltered_bmp
  fi
fi

#need for hgtfile etc.
frame=`pwd | rev | cut -d '/' -f3 | rev`
tr=`echo $frame | cut -c -3 | sed 's/^0//' | sed 's/^0//'`
hgtfile=$LiCSAR_public/$tr/$frame/metadata/$frame.geo.hgt.tif
hillshadefile=$LiCSAR_public/$tr/$frame/metadata/$frame.geo.hillshade.nc

if [ ! -f $hgtfile ]; then
  mkdir -p ../../geo
  hgtfile=../../geo/$frame.geo.hgt.tif
  if [ ! -f $hgtfile ]; then
   #try for the last time, download it..
   wget --no-check-certificate -O $hgtfile https://gws-access.ceda.ac.uk/public/nceo_geohazards/LiCSAR_products/$tr/$frame/metadata/$frame.geo.hgt.tif 2>/dev/null
  fi
  if [ ! -f $hgtfile ]; then
   echo "warning, no hgt tiff file found. there will be no hillshade"
  fi
fi

if [ ! -f $hillshadefile ]; then
 if [ ! -d $LiCSAR_public/$tr/$frame/metadata ]; then
  #in such case do it only locally
  hillshadefile='hillshade.nc'
 fi
fi

if [ $do_unw -eq 1 ]; then
  echo "picturing unwrapped ifg"
  unw_bmp=$pair.geo.unw.full.png
  rm $unw_bmp 2>/dev/null
  #finally let's get the unwrapped png:
  #so let's get 2% away.. or only 0.1...
  #minmaxcolour=`gmt grdinfo -T+a2 $unw_tif`
  #trimming 0.1% of all values--- this seems pretty enough!
  gmt grdclip $unw_tif -Gtemp.nc -Sr0/NaN
  #median=`gmt grdinfo -C -L1 temp.nc | gawk {'print $12'}`
  gmt grdmath temp.nc temp.nc MEDIAN SUB = temp3.nc
  mv temp3.nc temp.nc
  #minmaxcolour=`gmt grdinfo -T+a0.1+s temp.nc`
  minmaxcolour=`gmt grdinfo -T+a1+s temp.nc`
  #gmt makecpt -C$LiCSARpath/misc/colourmap.cpt $minmaxcolour/0.025 >unw.cpt 
  gmt makecpt -C$LiCSARpath/misc/colourmap.cpt -Iz $minmaxcolour/0.025 >unw.cpt
  
  hillshadecmd=''
  if [ ! -f $hillshadefile ]; then
   #do hillshade to unw..
   echo "(also generating hillshade)"
   
   if [ -f $hgtfile ]; then
    #gmt grdhisteq $hgtfile -Ghgt.grd -N
    #for ascending: 258 deg, descending: 102 deg
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
    hillshadecmd="-I"$hillshadefile
  fi
  #generate figure, incl. hillshade
  gmt grdimage temp.nc -Cunw.cpt -JM1 -Q -nn+t0.1 $hillshadecmd -A$unw_bmp 2>/dev/null
  if [ ! -f $unw_bmp ]; then
    echo "error with hillshade, trying to regenerate it"
    #gmt grdsample $hillshadefile -Gtemphill.nc -Rtemp.nc #2>/dev/null
    opass=`echo $frame | cut -c 4`
    if [ $opass == 'A' ]; then
      deg=258
    else
      deg=102
    fi
    gmt grdgradient -A$deg -Nt1 -Gtemphill.nc -Rtemp.nc $hgtfile
    gmt grdimage temp.nc -Cunw.cpt -JM1 -Q -nn+t0.1 -Itemphill.nc -A$unw_bmp
  fi
  if [ ! -f $unw_bmp ]; then
    echo "some error, perhaps with hillshade.. skipping it"
    gmt grdimage temp.nc -Cunw.cpt -JM1 -Q -nn+t0.1 -Atempunwbmp.png
    #in such case we should be ok with 255 colour palette
    convert tempunwbmp.png PNG8:$unw_bmp
    rm tempunwbmp.png
  fi
  
  #need to prepare a colorbar based on these values!!!!
  #you know... the rounding here is not really important... just a colourbar.. or not?
  mincol=`echo $minmaxcolour | cut -d '/' -f1 | cut -d 'T' -f2`
  maxcol=`echo $minmaxcolour | cut -d '/' -f2 `
  #expecting sentinel
  minval=`python -c 'print(round('$mincol'*5.546/(4*3.14159265)))'`
  maxval=`python -c 'print(round('$maxcol'*5.546/(4*3.14159265)))'`
  #add also real min and max values
  minmaxreal=`gmt grdinfo -T temp.nc`
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
  rm unw.cpt temp.nc temp_scale_unw.png hillsh.ti* 2>/dev/null

  
   #gmt grdimage pokus.tif -C$LiCSARpath/misc/pha.cpt -JM1 -nn+t0.1 -I hillsh.tif -A$unw_bmp
   #convert output3.png 20200907_20200913.geo.cc.full.png -compose CopyOpacity -composite PNG32:fixed.png
   #convert -flatten fixed.png fixed2.png
fi






rm gmt.history 2>/dev/null
#generating previews
cp -r $LiCSARpath/misc/kmlfiles files
#convert $unw_bmp files/`basename $unw_bmp .bmp`.png
#convert $ifg_bmp files/`basename $ifg_bmp .bmp`.png
cp $unw_bmp $ifg_bmp $coh_bmp files/.
#getting templates and fit it
cp $LiCSARpath/misc/template.kml $pair.kml



# update 2020: include unfiltered if exists
if [ -f $ifg_unfiltered_bmp ]; then
 #echo "including also unfiltered file"
 cp $ifg_unfiltered_bmp files/.
 cp $LiCSARpath/misc/template_unfilt.kml $pair.kml
fi

if [ -f scalebar_unwrapped.png ]; then
 cp scalebar_unwrapped.png files/.
fi


#sed -i 's///' $pair.kml
#sed -i 's/TYEART/'`date +'%Y'`'/' $pair.kml
sed -i 's/TFIRSTDATET/'`echo $pair | cut -d '_' -f1`'/' $pair.kml
sed -i 's/TLASTDATET/'`echo $pair | cut -d '_' -f2`'/' $pair.kml
#sed -i 's/TFILEUNWT/'`basename $unw_bmp .bmp`.png'/' $pair.kml
#sed -i 's/TFILEIFGT/'`basename $ifg_bmp .bmp`.png'/' $pair.kml
sed -i 's/TFILEUNWT/'`basename $unw_bmp`'/' $pair.kml
sed -i 's/TFILEIFGT/'`basename $ifg_bmp`'/' $pair.kml
sed -i 's/TFILECOHT/'`basename $coh_bmp`'/' $pair.kml

if [ -f $ifg_unfiltered_bmp ]; then
 sed -i 's/TFILEUNFIFGT/'`basename $ifg_unfiltered_bmp`'/' $pair.kml
fi

#coordinates
gdalinfo -stats $unw_tif | grep Coordinates -A4 | tail -n+2 | cut -d '(' -f2 | cut -d ')' -f1 > tmp.coord
rm tmp.ewsn 2>/dev/null
for x in `cat tmp.coord | cut -d ',' -f1 | sort -n -u`; do echo $x >> tmp.ewsn; done
for x in `cat tmp.coord | cut -d ',' -f2 | sort -n -u`; do echo $x >> tmp.ewsn; done

E=`sed '1q;d' tmp.ewsn`
W=`sed '2q;d' tmp.ewsn`
S=`sed '3q;d' tmp.ewsn`
N=`sed '4q;d' tmp.ewsn`
sed -i 's/TEASTT/'$E'/' $pair.kml
sed -i 's/TWESTT/'$W'/' $pair.kml
sed -i 's/TSOUTHT/'$S'/' $pair.kml
sed -i 's/TNORTHT/'$N'/' $pair.kml

centerlat=`python -c "print(("$S"+"$N")/2)"`
centerlon=`python -c "print(("$E"+"$W")/2)"`
sed -i 's/TCENLATT/'$centerlat'/' $pair.kml
sed -i 's/TCENLONT/'$centerlon'/' $pair.kml


#landmask
  if [ -f $LiCSAR_public/$tr/$frame/metadata/$frame.geo.landmask.tif ]; then
   maskfile=$LiCSAR_public/$tr/$frame/metadata/$frame.geo.landmask.tif
  elif [ -f $frame.geo.landmask.tif ]; then
   maskfile=$frame.geo.landmask.tif
  else 
   wget --no-check-certificate -O $frame.geo.landmask.tif https://gws-access.ceda.ac.uk/public/nceo_geohazards/LiCSAR_products/$tr/$frame/metadata/$frame.geo.landmask.tif 2>/dev/null
   if [ -f $frame.geo.landmask.tif ]; then maskfile=$frame.geo.landmask.tif; else maskfile=''; fi
  fi

if [ ! -z $maskfile ]; then
 maskmin=`gdalinfo -stats $maskfile 2>/dev/null | grep ICS_MINIMUM | cut -d '=' -f2`
 if [ $maskmin -eq 0 ]; then
  echo "applying landmask"
  #maskfilepng=`echo $maskfile | rev | cut -c 4- | rev`png
  maskfilepng=landmask.png
  if [ ! -f $maskfilepng ]; then
   gmt grdimage $maskfile -JM1 -R$pair.geo.diff_pha.tif -Cwhite,black -Gwhite -A$maskfilepng
  fi
  cp $maskfilepng tempmask.png
  
  for bmp in `ls files/$pair*png`; do
   bmpmasked=$bmp.masked.png;
   dims=`identify $bmp 2>/dev/null | gawk {'print $3'}`
   if [ ! `identify tempmask.png | gawk {'print $3'}` == $dims ]; then
    convert $maskfilepng -resize $dims\! tempmask.png
   fi
   #from  https://stackoverflow.com/questions/59935181/mask-png-image-with-black-and-white-mask
   #this is all what is needed for unw pictures
   if [ $bmp == files/$unw_bmp ]; then
    convert $bmp \( tempmask.png -negate \) \( -clone 0 -transparent black +transparent black \) -insert 0 -composite tempp.png
    convert tempp.png PNG8:$bmpmasked
    rm tempp.png
   else
    if [ ! -f tempcoh.png ]; then
     if [ ! `identify $coh_bmp | gawk {'print $3'}` == $dims ]; then
      convert $coh_bmp -resize $dims\! tempcoh.png
     else
      cp $coh_bmp tempcoh.png
     fi
    fi
    convert $bmp tempcoh.png -compose CopyOpacity -composite PNG32:fixed.png
    #convert -flatten fixed.png fixed2.png
    convert fixed.png \( tempmask.png -negate \) \( -clone 0 -transparent black +transparent black \) -insert 0 -composite tempp.png
   #convert -flatten fixed.png fixed2.png
    convert -flatten tempp.png tempp2.png
    convert -transparent white tempp2.png PNG8:$bmpmasked
    rm tempp.png tempp2.png fixed.png tempcoh.png;
   fi
   mv $bmpmasked $bmp;
  done
  rm tempmask.png #$maskfilepng
 fi
fi

zip -r $pair.zip $pair.kml files >/dev/null
if [ $opt == 0 ]; then
 outfile=$pair.kmz
else
 outfile=$opt'_'$pair.kmz
fi
mv $pair.zip $outfile



#now this is to update the previews...
if [ $update_previews -eq 1 ]; then
 echo "updating previews"
 #cp $LiCSARpath/misc/kmlfiles/scalebar_wrapped.png $LiCSARpath/misc/kmlfiles/scalebar_coh.png .
 #cp scalebar_wrapped.png scalebar_wrapped2.png
 convert $unw_bmp -resize 680x \( scalebar_unwrapped.png -resize 385x  -background none -gravity center \) -gravity southwest -geometry +7+7 -composite -flatten -transparent black $pair.geo.unw.png
 #convert $ifg_bmp -resize 400x \( scalebar_wrapped.png -resize 385x  -background none -gravity center \) -gravity southwest -geometry +7+10 -composite -flatten -transparent black $pair.geo.diff.png
 convert $ifg_bmp -resize 680x -flatten -transparent black $pair.geo.diff.png
 if [ -f $ifg_unfiltered_bmp ]; then
  #convert $ifg_unfiltered_bmp -resize 400x \( scalebar_wrapped2.png -resize 385x -background none -gravity center \) -gravity southwest -geometry +7+10 -composite -flatten -transparent black $pair.geo.diff_unfiltered.png
  convert $ifg_unfiltered_bmp -resize 680x -flatten -transparent black $pair.geo.diff_unfiltered.png
 fi
 #convert $coh_bmp -resize 400x \( scalebar_coh.png -resize 385x  -background none -gravity center \) -gravity southwest -geometry +7+10 -composite -flatten PNG8:$pair.geo.cc.png
 convert $coh_bmp -resize 680x -flatten PNG8:$pair.geo.cc.png
 rm scalebar_* 
fi


rm tmp.coord tmp.ewsn
rm -r files $pair.kml *.aux.xml *.full.png
rm landmask.png gmt.history 2>/dev/null
#actually also delete the 'full' ifgs - why should we use them? ...
if [ `echo $1 | grep -c LiCSAR_products` -gt 0 ]; then
 echo "done. Being in public folder, the file should be accessible at:"
 echo "https://gws-access.ceda.ac.uk/public/nceo_geohazards/"`echo $1 | cut -c 43-`"/"$outfile
fi

