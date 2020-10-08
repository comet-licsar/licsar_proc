#!/bin/bash
#this will force regenerate full previews
overwrite=1

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
 echo "generating fullscale preview"
 #doing also unw for full scale then
 unw_bmp=$pair.geo.unw.full.png
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
 #to flatten it:
 convert $ifg_bmp temptemp.png
 mv temptemp.png $ifg_bmp
fi

if [ $do_ifgu -eq 1 ]; then
 #perhaps we are to add unfiltered one:
  if [ -f $pair.geo.diff_unfiltered_pha.tif ]; then
    echo "picturing unfiltered ifg"
    gmt grdimage $pair.geo.diff_unfiltered_pha.tif -C$LiCSARpath/misc/pha.cpt -Q -JM1 -nn+t0.1 -A$ifg_unfiltered_bmp
    convert $ifg_unfiltered_bmp temptemp.png
    mv temptemp.png $ifg_unfiltered_bmp
  fi
fi

#need for hgtfile etc.
frame=`pwd | rev | cut -d '/' -f3 | rev`
tr=`echo $frame | cut -c -3 | sed 's/^0//' | sed 's/^0//'`
hgtfile=$LiCSAR_public/$tr/$frame/metadata/$frame.geo.hgt.tif

if [ $do_unw -eq 1 ]; then
  echo "picturing unwrapped ifg"
  unw_bmp=$pair.geo.unw.full.png
  #finally let's get the unwrapped png:
  #so let's get 2% away..
  minmaxcolour=`gmt grdinfo -T+a2 $unw_tif`
  gmt makecpt -C$LiCSARpath/misc/colourmap.cpt -Iz $minmaxcolour/0.025 >unw.cpt 
  #changing 0 to NaN - donno why, always problem with it
  gmt grdclip $unw_tif -Gtemp.nc -Sr0/NaN
  
  #do hillshade to unw..
  echo "(generating hillshade)"
  hillshadecmd=''
  if [ -f $hgtfile ]; then
    #gmt grdhisteq $hgtfile -Ghgt.grd -N
    #for ascending: 258 deg, descending: 102 deg
    opass=`echo $frame | cut -c 4`
    if [ $opass == 'A' ]; then
      deg=258
    else
      deg=102
    fi
    gmt grdgradient -A$deg -Nt1 -Ghillsh.tif=gd:GTiff $hgtfile
    hillshadecmd="-Ihillsh.tif"
  fi
  #generate figure, incl. hillshade
  gmt grdimage temp.nc -Cunw.cpt -JM1 -Q -nn+t0.1 $hillshadecmd -A$unw_bmp
  #need to prepare a colorbar based on these values!!!!
  #you know... the rounding here is not really important... just a colourbar.. or not?
  mincol=`echo $minmaxcolour | cut -d '/' -f1 | cut -d 'T' -f2`
  maxcol=`echo $minmaxcolour | cut -d '/' -f2 `
  #expecting sentinel
  minval=`python -c 'print(round('$mincol'*5.546/(4*3.14159265)))'`
  maxval=`python -c 'print(round('$maxcol'*5.546/(4*3.14159265)))'`
  #burn them to the scalebar
  convert -font helvetica -fill black -pointsize 40 -draw "text 15,115 '"$minval"'" $LiCSARpath/misc/scalebar_unwrapped_empty.png temp_scale_unw.png
  convert -font helvetica -fill black -pointsize 40 -draw "text 1100,115 '"$maxval" cm'" temp_scale_unw.png scalebar_unwrapped.png
  rm unw.cpt temp.nc temp_scale_unw.png hillsh.ti* 2>/dev/null

  
   #gmt grdimage pokus.tif -C$LiCSARpath/misc/pha.cpt -JM1 -nn+t0.1 -I hillsh.tif -A$unw_bmp
   #convert output3.png 20200907_20200913.geo.cc.full.png -compose CopyOpacity -composite PNG32:fixed.png
   #convert -flatten fixed.png fixed2.png
fi




#~ if [ ! -f $unw_bmp ]; then
#~ #  echo "generating fullscale preview"
   #~ gmt grdmath $unw_tif 3 DIV WRAP = pokus.tif=gd:GTiff
   #~ gmt grdimage pokus.tif -C$LiCSARpath/misc/pha.cpt -JM1 -nn+t0.1 -A$unw_bmp
   #~ #hgtshade
   #~ hgtfile=$LiCSAR_public/$tr/$frame/metadata/$frame.geo.hgt.tif
   #~ ...........
   #~ gmt grdhisteq $hgtfile -Ghgt.grd -N
   #~ #for ascending: 258 deg
   #~ gmt grdgradient -A258 -Nt1 -Ghillsh.tif=gd:GTiff $hgtfile #hgt.grd
   #~ gmt grdimage pokus.tif -C$LiCSARpath/misc/pha.cpt -JM1 -nn+t0.1 -I hillsh.tif -A$unw_bmp
   #~ convert output3.png 20200907_20200913.geo.cc.full.png -compose CopyOpacity -composite PNG32:fixed.png
   #~ convert -flatten fixed.png fixed2.png

   #~ #gmt grdhisteq hillsh.tif -GhillshN.tif=gd:GTiff -N
   #~ #gmt grdimage hillshN.tif -Cgray -JM1 -P -Ahillsh.png
#~ #grdimage $in_grd -Jx0.25c -O -K -C$in_cpt -I$in_shadow -Y3c >> $
   #~ rm pokus.tif
#~ #  gmt grdimage $pair.geo.diff_pha.tif -C$LiCSARpath/misc/unw.cpt -JM1 -A$unw_bmp
  #~ #convert -transparent black temp_blk.png $ifg_bmp
  #~ #rm temp_blk.png
#~ fi

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
 echo "including also unfiltered file"
 cp $ifg_unfiltered_bmp files/.
 cp $LiCSARpath/misc/template_unfilt.kml $pair.kml
fi

if [ -f scalebar_unwrapped.png ]; then
 mv scalebar_unwrapped.png files/.
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
 maskmin=`gdalinfo -stats $maskfile | grep ICS_MINIMUM | cut -d '=' -f2`
 if [ $maskmin -eq 0 ]; then
  echo "applying landmask"
  maskfilepng=`echo $maskfile | rev | cut -c 4- | rev`png
  if [ ! -f $maskfilepng ]; then
   gmt grdimage $maskfile -JM1 -Cwhite,black -Gwhite -A$maskfilepng
  fi
  cp $maskfilepng tempmask.png
  for bmp in `ls files/$pair*png`; do
   bmpmasked=$bmp.masked.png;
   dims=`identify $bmp 2>/dev/null | gawk {'print $3'}`
   if [ ! `identify tempmask.png | gawk {'print $3'}` == $dims ]; then
    convert $maskfilepng -resize $dims\! tempmask.png
   fi
   #from  https://stackoverflow.com/questions/59935181/mask-png-image-with-black-and-white-mask
   convert $bmp \( tempmask.png -negate \) \( -clone 0 -transparent black +transparent black \) -insert 0 -composite tempp.png
   #this is all what is needed for unw pictures
   if [ $bmp == files/$unw_bmp ]; then
    mv tempp.png $bmpmasked
   else
    convert -flatten tempp.png tempp2.png
    convert -transparent white tempp2.png $bmpmasked
    rm tempp.png tempp2.png;
   fi
   mv $bmpmasked $bmp;
  done
  rm tempmask.png
 fi
fi

zip -r $pair.zip $pair.kml files >/dev/null
if [ $opt == 0 ]; then
 outfile=$pair.kmz
else
 outfile=$opt'_'$pair.kmz
fi
mv $pair.zip $outfile

rm tmp.coord tmp.ewsn
rm -r files $pair.kml $unw_tif.aux.xml
if [ `echo $1 | grep -c LiCSAR_products` -gt 0 ]; then
 echo "done. Being in public folder, the file should be accessible at:"
 echo "https://gws-access.ceda.ac.uk/public/nceo_geohazards/"`echo $1 | cut -c 43-`"/"$outfile
fi

