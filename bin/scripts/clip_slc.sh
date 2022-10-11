#!/bin/bash
#this will clip the SLCs ...
# you may get the lonlats as:
# hgt=$LiCSAR_public/`track_from_frame $frame`/$frame/metadata/$frame.geo.hgt.tif
# gdallocationinfo $hgt $lon $lat
# so... if i have lon1, lon2 etc., i can just do... nah, let's do it through python, i.e.:
# import xarray as xr; import rioxarray
# a=rioxarray.open_rasterio(hgt)
# lon=-71.377; lat=-36.863; radius_km=25/2; radius_deg=radius_km/111
# a.sel(lon=(lon-radius_deg, lon+radius_deg), lat=(lat+radius_deg, lat-radius_deg))
# medhgt=float(a.sel(x=(lon-radius_deg, lon+radius_deg), y=(lat+radius_deg, lat-radius_deg), method='nearest').median())
# print(str(lon-radius_deg), lon+radius_deg, lat-radius_deg, lat+radius_deg, medhgt)

if [ -z $7 ]; then echo "parameters are:";
echo "clip_slc.sh SLCFOLDER OUTFOLDER lon1 lon2 lat1 lat2 hei"
echo "so e.g. clip_slc.sh RSLC/20200112 CLIPPED -28.36 -27.3 38.49 38.8 600"
echo "e.g. /gws/nopw/j04/nceo_geohazards_vol1/public/LiCSAR_products/58/058A_05279_131311/interferograms/20190317_20190323"
echo "second (optional) parameter is a prefix, instead of frame name"
echo "third (optional) parameter is an overwrite control - by default, this is 0 (use 1 to overwrite existing kmz)"
exit;
fi

slcpath=$1
epoch=`echo $slcpath | rev | cut -d '/' -f1 | rev`
slc=`ls $slcpath/$epoch.slc $slcpath/$epoch.rslc 2>/dev/null`
slcpar=$slc.par
if [ ! -f $slcpar ]; then echo "the folder "$1" seems empty, or no mosaic exists - exiting"; exit; fi


outdir=$2
mkdir -p $outdir 2>/dev/null

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
