#!/bin/bash
#this will clip the SLCs ...


if [ -z $7 ]; then echo "parameters are:";
echo "clip_slc.sh SLCFOLDER OUTFOLDER lon1 lon2 lat1 lat2 hei"
echo "so e.g. clip_slc.sh RSLC/20200112 CLIPPED -28.36 -27.3 38.49 38.8 600"
echo "e.g. /gws/nopw/j04/nceo_geohazards_vol1/public/LiCSAR_products/58/058A_05279_131311/interferograms/20190317_20190323"
echo "second (optional) parameter is a prefix, instead of frame name"
echo "third (optional) parameter is an overwrite control - by default, this is 0 (use 1 to overwrite existing kmz)"
exit;
fi

epoch=`echo $1 | rev | cut -d '/' -f1 | rev`
slcpath=$1
slc=$1/$epoch.slc $1/$epoch.rslc 2>/dev/null
if [ -z $slcpar ]; then echo "the folder "$1" seems empty, or no mosaic exists - exiting"; exit; fi
slcpar=$slc.par

outdir=$2
mkdir -p $2 2>/dev/null

lon1=$3
lon2=$4
lat1=$5
lat2=$6
hei=$7

coord_to_sarpix $slcpar - - $lat1 $lon1 $hei | grep "SLC/MLI range, azimuth pixel (int)" > corners_clip.tmp
coord_to_sarpix $slcpar - - $lat2 $lon2 $hei | grep "SLC/MLI range, azimuth pixel (int)" >> corners_clip.tmp
coord_to_sarpix $slcpar - - $lat1 $lon2 $hei | grep "SLC/MLI range, azimuth pixel (int)" >> corners_clip.tmp
coord_to_sarpix $slcpar - - $lat2 $lon1 $hei | grep "SLC/MLI range, azimuth pixel (int)" >> corners_clip.tmp

azi1=`cat corners_clip.tmp | rev | gawk {'print $1'} | rev | sort -n | head -n1`
azi2=`cat corners_clip.tmp | rev | gawk {'print $1'} | rev | sort -n | tail -n1`
let azidiff=azi2-azi1+1

rg1=`cat corners_clip.tmp | rev | gawk {'print $2'} | rev | sort -n | head -n1`
rg2=`cat corners_clip.tmp | rev | gawk {'print $2'} | rev | sort -n | tail -n1`
let rgdiff=rg2-rg1+1

# ok, now clip the mosaic
mkdir -p $outdir/$slcpath
echo "clipping "$slcpath
SLC_copy $slc $slcpar $outdir/$slc $outdir/$slcpar - - $rg1 $rgdiff $azi1 $azidiff - - >/dev/null
