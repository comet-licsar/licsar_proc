#!/bin/bash
#this will force regenerate full previews
overwrite=0
update_previews=0

if [ $overwrite -eq 1 ]; then
 update_previews=1
fi

if [ -z $1 ]; then echo "parameter is full path to the folder containing geotiffs";
echo "e.g. /gws/nopw/j04/nceo_geohazards_vol1/public/LiCSAR_products/58/058A_05279_131311/interferograms/20190317_20190323"
echo "second (optional) parameter is a prefix, instead of frame name"
echo "third (optional) parameter is an overwrite control - by default, this is 0 (use 1 to overwrite existing kmz)"
exit;
fi

if [ ! -z $2 ]; then opt=$2; else opt=0; fi
if [ ! -z $3 ]; then overwrite=$3; fi

pair=`basename $1`
if [ ! `echo $pair | wc -m` -eq 18 ]; then 
 echo "wrong path - it should end on pairname, e.g. 20190101_20190107"
 exit
fi

if [ ! -f $1/$pair.geo.diff_pha.tif ] || [ ! -f $1/$pair.geo.unw.tif ]; then
 echo "some of the geocoded product does not exist here"
 exit
fi







### let's make it
echo "version 10/2020 by ML - use of GMT GRD2KML"
echo "generating kmz file for pair "$pair

cd $1





frame=`pwd | rev | cut -d '/' -f3 | rev`
tr=`echo $frame | cut -c -3 | sed 's/^0//' | sed 's/^0//'`
echo "(frame "$frame " )"

if [ $opt == 0 ]; then
 outfile=$frame'_'$pair
 opt=$frame
else
 outfile=$opt'_'$pair
fi

if [ $overwrite -eq 0 ]; then
 if [ -f $outfile ]; then
  echo "the kmz already exists, not overwriting (change it by adding extra parameter)"
  exit
 fi
fi

unw_tif=$pair.geo.unw.tif
if [ -f $pair.geo.diff_unfiltered_pha.tif ]; then
 do_ifgu=1
else
 do_ifgu=0
fi

#need for hillshade
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

if [ ! -f $hillshadefile ]; then
   #do hillshade to unw..
   echo "(generating hillshade)"
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





echo "applying landmask if needed"
#check and/or regenerate landmask
if [ -f $LiCSAR_public/$tr/$frame/metadata/$frame.geo.landmask.tif ]; then
   maskfile=$LiCSAR_public/$tr/$frame/metadata/$frame.geo.landmask.tif
elif [ -f $frame.geo.landmask.tif ]; then
   maskfile=$frame.geo.landmask.tif
else 
   wget --no-check-certificate -O $frame.geo.landmask.tif https://gws-access.ceda.ac.uk/public/nceo_geohazards/LiCSAR_products/$tr/$frame/metadata/$frame.geo.landmask.tif 2>/dev/null
   if [ -f $frame.geo.landmask.tif ]; then maskfile=$frame.geo.landmask.tif; else maskfile=''; fi
fi

landmask=0
if [ ! -z $maskfile ]; then
 maskmin=`gdalinfo -stats $maskfile 2>/dev/null | grep ICS_MINIMUM | cut -d '=' -f2`
 if [ $maskmin -eq 0 ]; then
  #resample to the AOI (fix for the missing bursts etc.)
  gmt grdsample $maskfile -Gtempmask.nc -R$unw_tif -nn+t0.1 2>/dev/null
  #gmt grdcut -N -R$unw_tif -Gtempmask.nc $maskfile 2>/dev/null
  if [ -f tempmask.nc ]; then
   echo "will apply landmask"
   masknc=tempmask.nc 
   landmask=1
  fi
 fi
fi

#we will do a parsed kml file
kmlfile=doc.kml
center=`gdalinfo $unw_tif | grep Center`
centerlon=`echo $center | gawk {'print $3'} | cut -d ',' -f1`
centerlat=`echo $center | gawk {'print $4'} | cut -d ')' -f1`
cat << EOF > $kmlfile
<?xml version="1.0" encoding="UTF-8"?>
<kml xmlns="http://www.opengis.net/kml/2.2" xmlns:gx="http://www.google.com/kml/ext/2.2" xmlns:kml="http://www.opengis.net/kml/2.2" xmlns:atom="http://www.w3.org/2005/Atom">
<Folder>
        <name>$frame : $pair</name>
        <description>$pair</description>
        <open>1</open>
          <LookAt>
           <longitude>$centerlon</longitude>
           <latitude>$centerlat</latitude>
           <altitude>0</altitude>
           <heading>0</heading>
           <tilt>0</tilt>
           <range>500000</range>
         </LookAt>
    <Folder>
        <Style id="radioLayer">
         <ListStyle>
           <listItemType>radioFolder</listItemType>
         </ListStyle>
        </Style>
        <name>InSAR data</name>
            <open>1</open>
EOF


#the sorting here matters!
todo="diff_pha cc unw"
tododirs=''
if [ $do_ifgu -eq 1 ]; then
 todo="diff_pha diff_unfiltered_pha cc unw"
fi

for topic in $todo; do
 if [ $topic == "cc" ]; then
  keyword="coherence"
  cpt=$LiCSARpath/misc/cc.cpt
  scalebarfile=$LiCSARpath/misc/kmlfiles/scalebar_coh.png
  extracmd=''
  tododirs=$tododirs" "$keyword
  visible=0
 elif [ $topic == "diff_pha" ]; then
  keyword="filtered_ifg"
  cpt=$LiCSARpath/misc/pha.cpt
  scalebarfile=$LiCSARpath/misc/kmlfiles/scalebar_wrapped.png
  extracmd=''
  tododirs=$tododirs" "$keyword
  visible=0
 elif [ $topic == "diff_unfiltered_pha" ]; then
  keyword="unfiltered_ifg"
  cpt=$LiCSARpath/misc/pha.cpt
  scalebarfile=$LiCSARpath/misc/kmlfiles/scalebar_wrapped.png
  extracmd=''
  tododirs=$tododirs" "$keyword
  visible=0
 elif [ $topic == "unw" ]; then
  keyword="unwrapped_ifg"
  #i needed to do hillshade resample later on... so only:
  extracmd=''
#  extracmd=$hillshadecmd
  #have to generate the unw cpt..
  cpt=''
  tododirs=$tododirs" "$keyword
  visible=1
 fi
 
 echo "picturing "$keyword
 infile=$pair.geo.$topic.tif
 if [ $landmask -eq 1 ]; then
  #should apply landmask  -  file $maskfile
  rm tempmasked.nc 2>/dev/null
  gmt grdclip $infile -Gtempinfile1.nc -SrNaN/0
  gmt grdsample tempinfile1.nc -Gtempinfile2.nc -R$masknc -nn+t0.1 2>/dev/null
  gmt grdmath tempinfile2.nc $masknc MUL = temp2.nc
  gmt grdclip temp2.nc -Gtempmasked.nc -Sr0/NaN
  rm temp2.nc 2>/dev/null
  if [ ! -f tempmasked.nc ]; then
   echo "problem with landmasking, keeping without mask"
  else
   infile=tempmasked.nc
  fi
 fi
 
 #sometimes this does not work.. so we need to just.. change 0 for NaNs..
 if [ ! $infile == "tempmasked.nc" ]; then
  rm temp.nc 2>/dev/null
  gmt grdclip $infile -Gtemp.nc -Sr0/NaN
  if [ -f temp.nc ]; then
   infile=temp.nc
  fi
 fi
 

 if [ $topic == "unw" ]; then
  #for unwrapped data, i need to fix by the median:
  gmt grdmath $infile $infile MEDIAN SUB = temp3.nc
  mv temp3.nc temp.nc
  if [ -f temp.nc ]; then
   infile=temp.nc
  fi
  
  #resample hillshade for unw data only
  hillshadecmd=''
  if [ -f $hillshadefile ]; then
      hillshadenc='temphillshade.nc'
      rm $hillshadenc 2>/dev/null
      gmt grdsample $hillshadefile -G$hillshadenc -R$infile -nn+t0.1 2>/dev/null
  fi
  if [ -f $hillshadenc ]; then
      hillshadecmd="-I"$hillshadenc
      extracmd=$hillshadecmd
  else
      echo "error with hillshade - skipping it"
      extracmd=''
  fi
  #create cpt - trim 1% off from both sides of histogram
  minmaxcolour=`gmt grdinfo -T+a1+s $infile`
  gmt makecpt -C$LiCSARpath/misc/colourmap.cpt -Iz $minmaxcolour/0.025 >unw.cpt
  cpt='unw.cpt'
  
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
  
  scalebarfile='scalebar_unwrapped.png'
  
 fi
 
 #now generating kml from it
 gmt grd2kml -Ag -C$cpt -nn+t0.1 -T$keyword -N$keyword $extracmd $infile 2>/dev/null
 if [ ! -f $keyword/$keyword.kml ]; then
  echo "some error - perhaps with hillshade.."
  rm -r $keyword
  gmt grd2kml -Ag -C$cpt -nn+t0.1 -T$keyword -N$keyword $infile 2>/dev/null
 fi
 #compressing to PNG8
 for x in `ls $keyword/*/*png`; do mv $x tt.png; convert tt.png PNG8:$x; done; rm tt.png

 
 #now edit the kml file - add scalebar there
 cp $scalebarfile $keyword/scalebar.png
 
cat << EOF >> $kmlfile
<Document>
<name>$keyword</name>
<Document>
    <name>$keyword</name>
<Style>
EOF
 
 grep "<Style>" -A 99999999999 $keyword/$keyword.kml | tail -n+2 | head -n-2 | sed 's/<href>/<href>'$keyword'\//' >> $kmlfile
 
 
cat << EOF >> $kmlfile
 </Document>
 <Document>
    <name>colour bar</name>
    <ScreenOverlay>
        <name>colour bar</name>
        <visibility>1</visibility>
        <Icon>
          <href>$keyword/scalebar.png</href>
        </Icon>
        <overlayXY x="0" y="0" xunits="fraction" yunits="fraction"/>
        <screenXY x="0" y="0" xunits="fraction" yunits="fraction"/>
        <rotationXY x="0" y="0" xunits="fraction" yunits="fraction"/>
        <size x="500" y="0" xunits="pixels" yunits="fraction"/>
    </ScreenOverlay>
 </Document>
</Document>
EOF

done



#finalise the KML
echo "adding logos"
mkdir logos
cp $LiCSARpath/misc/kmlfiles/lics.png logos/.
cp $LiCSARpath/misc/kmlfiles/nerc.png logos/.
cp $LiCSARpath/misc/kmlfiles/comet.png logos/.

cat << EOF >> $kmlfile
</Folder>

<Folder>
<name>Logos</name>
<ScreenOverlay>
        <name>LiCS NERC</name>
        <Icon>
          <href>logos/lics.png</href>
        </Icon>
        <overlayXY x="0" y="1" xunits="fraction" yunits="fraction"/>
        <screenXY x="0" y="1" xunits="fraction" yunits="fraction"/>
        <rotationXY x="0" y="0" xunits="fraction" yunits="fraction"/>
        <size x="0.12" y="0" xunits="fraction" yunits="fraction"/>
</ScreenOverlay>
<ScreenOverlay>
        <name>COMET NERC</name>
        <Icon>
          <href>logos/comet.png</href>
        </Icon>
        <overlayXY x="0.9" y="1" xunits="fraction" yunits="fraction"/>
        <screenXY x="0.9" y="1" xunits="fraction" yunits="fraction"/>
        <rotationXY x="0" y="0" xunits="fraction" yunits="fraction"/>
        <size x="0.16" y="0" xunits="fraction" yunits="fraction"/>
</ScreenOverlay>
<ScreenOverlay>
        <name>NERC</name>
        <Icon>
          <href>logos/nerc.png</href>
        </Icon>
        <overlayXY x="0.8" y="0.0" xunits="fraction" yunits="fraction"/>
        <screenXY x="0.8" y="0.025" xunits="fraction" yunits="fraction"/>
        <rotationXY x="0" y="0" xunits="fraction" yunits="fraction"/>
        <size x="0.085" y="0" xunits="fraction" yunits="fraction"/>
</ScreenOverlay>
</Folder>

</Folder>
</kml>
EOF

7za a $outfile.zip logos $tododirs $kmlfile >/dev/null 2>/dev/null
mv $outfile.zip $outfile.kmz


#cleaning
rm -r $tododirs logos gmt.history temp*nc temp_scale_unw.png scalebar_unwrapped.png unw.cpt $kmlfile 2>/dev/null
