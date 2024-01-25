#!/bin/bash
#this will force regenerate full previews
overwrite=0
update_previews=0
source $LiCSARpath/lib/LiCSAR_bash_lib.sh


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
#echo "version 10/2020 by ML - use of GMT GRD2KML"
echo "version 02/2023 by ML - adding rng/azi offsets if exist"
echo "generating kmz file for pair "$pair

cd $1





frame=`pwd | rev | cut -d '/' -f3 | rev`
tr=`echo $frame | cut -c -3 | sed 's/^0//' | sed 's/^0//'`
echo "(frame "$frame " )"

if [ $opt == 0 ]; then
 outfilee=$frame'_'$pair
 opt=$frame
else
 outfilee=$opt'_'$pair
fi

if [ $overwrite -eq 0 ]; then
 if [ -f $outfilee.kmz ]; then
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

if [ -f $pair.geo.rng.tif ]; then
 do_off=1
 echo "including offsets - not adding unfiltered ifgs"
 do_ifgu=0
else
 do_off=0
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
EOF

# to add USGS eq event kml - needs to be in the folder of the ifg
#eqkml=`ls -t ../../*.kml 2>/dev/null | head -n1 2>/dev/null`
eqkml=`ls -t *.kml 2>/dev/null | head -n1 2>/dev/null`
if [ ! -z $eqkml ]; then
 if [ `tail -n1 $eqkml | grep -c USGS` -gt 0 ]; then
  tail -n1 $eqkml | cut -d '>' -f 2- | rev | cut -d '<' -f 2- | rev >> $kmlfile
 fi
fi

cat << EOF >> $kmlfile
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
if [ $do_off -eq 1 ]; then
 todo="diff_pha cc unw azi rng"
fi

for topic in $todo; do
 if [ $topic == "cc" ]; then
  keyword="coherence"
  cpt=$LiCSARpath/misc/cc.cpt
  scalebarfile=$LiCSARpath/misc/kmlfiles/scalebar_coh.png
  extracmd=''
  tododirs=$tododirs" "$keyword
  visible=0
  maskit=0
 elif [ $topic == "diff_pha" ]; then
  keyword="filtered_ifg"
  cpt=$LiCSARpath/misc/pha.cpt
  scalebarfile=$LiCSARpath/misc/kmlfiles/scalebar_wrapped.png
  extracmd=''
  tododirs=$tododirs" "$keyword
  visible=0
  maskit=1
 elif [ $topic == "diff_unfiltered_pha" ]; then
  keyword="unfiltered_ifg"
  cpt=$LiCSARpath/misc/pha.cpt
  scalebarfile=$LiCSARpath/misc/kmlfiles/scalebar_wrapped.png
  extracmd=''
  tododirs=$tododirs" "$keyword
  visible=0
  maskit=0
 elif [ $topic == "unw" ]; then
  keyword="unwrapped_ifg"
  create_preview_unwrapped $unw_tif $frame 1
  #cpt=unw.cpt
  scalebarfile=scalebar_unwrapped.png
  tododirs=$tododirs" "$keyword
  visible=1
  #maskit is 0 here because it was already done by create_preview_unwrapped...
  maskit=0
 elif [ $topic == "rng" ]; then
  keyword="range_offsets"
  create_preview_offsets `ls $pair.geo.rng.tif` $frame 10 1
  cpt=$topic'off'.cpt
  scalebarfile=scalebar_rng.png
  tododirs=$tododirs" "$keyword
  visible=0
  #maskit is 0 here because it was already done...
  maskit=0
 elif [ $topic == "azi" ]; then
  keyword="azimuth_offsets"
  create_preview_offsets `ls $pair.geo.azi.tif` $frame 10 1
  cpt=$topic'off'.cpt
  scalebarfile=scalebar_azi.png
  tododirs=$tododirs" "$keyword
  visible=0
  #maskit is 0 here because it was already done...
  maskit=0
 fi
 
 infile=$pair.geo.$topic.tif
 #now generating kml from it
 
 #gmt grd2kml -Ag -C$cpt -nn+t0.1 -T$keyword -N$keyword $extracmd $infile 2>/dev/null
 if [ ! -f $keyword/$keyword.kml ]; then
  #echo "some error - perhaps with hillshade.."
  echo "preparing kml layer for "$keyword
  rm -r $keyword 2>/dev/null
  if [ $maskit -eq 1 ]; then 
   maskedf=`prepare_landmask $infile $frame`
   if [ ! -z $maskedf ]; then infile=$maskedf; fi
  fi
  #this is necessary for tif, but not for maskedf.. but why not
  gmt grdclip $infile -Gtokml.nc -Sr0/NaN
  gmt grd2kml -Ag -C$cpt -nn+t0.1 -T$keyword -N$keyword tokml.nc 2>/dev/null
  rm tokml.nc $maskedf 2>/dev/null
 fi
 #compressing to PNG8
 #for x in `ls $keyword/*/*png`; do mv $x tt.png; convert tt.png -transparent black PNG8:$x; done; rm tt.png
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

7za a $outfilee.zip logos $tododirs $kmlfile >/dev/null 2>/dev/null
mv $outfilee.zip $outfilee.kmz


#cleaning
rm -r $tododirs logos gmt.history scalebar_*.png *.cpt $kmlfile 2>/dev/null
