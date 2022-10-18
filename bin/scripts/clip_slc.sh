#!/bin/bash
#this will clip the SLCs ...

# you may get hgt value as:
# hgt=$LiCSAR_public/`track_from_frame $frame`/$frame/metadata/$frame.geo.hgt.tif
# gdallocationinfo $hgt $lon $lat

# but better do it through python, i.e.:
# import xarray as xr; import rioxarray; import os
# frame='083D_12636_131313'
# hgt=os.path.join(os.environ['LiCSAR_public'], str(int(frame[:3])), frame, 'metadata', frame+'.geo.hgt.tif')
# a=rioxarray.open_rasterio(hgt)
# lon=-71.377; lat=-36.863; radius_km=25/2; radius_deg=radius_km/111
# # a.sel(lon=(lon-radius_deg, lon+radius_deg), lat=(lat+radius_deg, lat-radius_deg))
# medhgt=float(a.sel(x=(lon-radius_deg, lon+radius_deg), y=(lat+radius_deg, lat-radius_deg), method='nearest').median())
# print(str(lon-radius_deg), lon+radius_deg, lat-radius_deg, lat+radius_deg, medhgt)

# now we can clip all RSLCs (if they do not exist)
# for x in `ls RSLC/* -d`; do 
# if [ ! -d $x ]; then clip_slc.sh $x $outdir $lon1 $lon2 $lat1 $lat2 $hei; fi
# done

# and then generate hires geo (if it doesn't exist)
# cd volclip
# DEMDIR=$LiCSAR_proc..../DEM
# geodir='geo'
# masterslcdir='RSLC/'$master
# outres=0.00027027 # for 15 m, or 30/111000 for 30 m res.
# in python:
# python3 -c "from LiCSAR_lib.coreg_lib import geocode_dem; \
  geocode_dem('"$masterslcdir"', '"$geodir"', '"$DEMDIR"' , '.', '"$master"', "$outres")"

# ok, now time to generate ifgs and unws

if [ -z $7 ]; then echo "parameters are:";
echo "clip_slc.sh OUTFOLDER lon1 lon2 lat1 lat2 hei resolution"
echo "so e.g. clip_slc.sh CLIPPED -28.36 -27.3 38.49 38.8 600 0.00027"
exit;
fi

if [ -d RSLC ]; then echo "you need to be in the frame proc folder, i.e. the one with RSLC folder"; exit; fi

master=`basename RSLC/../geo/20??????.hgt | cut -d '.' -f1`
slc=`ls RSLC/$master/$master.rslc`
slcpar=$slc.par
if [ ! -f $slcpar ]; then echo "the folder "$1" seems empty, or no mosaic exists - exiting"; exit; fi


outdir=$1
mkdir -p $outdir 2>/dev/null

lon1=$2
lon2=$3
lat1=$4
lat2=$5
hei=$6
resol=$7
rgl=`echo $resol"*111000/2.3" | bc`
azl=`echo $resol"*111000/14" | bc`

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

# ok, now clip the mosaics

mkdir -p $outdir/RSLC $outdir/log 2>/dev/null
for x in `ls RSLC | grep 20`; do 
 if [ ! -d $outdir/RSLC/$x ]; then
   echo "clipping "$x
   mkdir -p $outdir/RSLC/$x
   SLC_copy RSLC/$x/$x.rslc RSLC/$x/$x.rslc.par $outdir/RSLC/$x/$x.rslc $outdir/RSLC/$x/$x.rslc.par - - $rg1 $rgdiff $azi1 $azidiff - - >/dev/null
   multi_look RSLC/$x/$x.rslc RSLC/$x/$x.rslc.par RSLC/$x/$x.mli RSLC/$x/$x.mli.par $rgl $azl
   # create_geoctiffs_to_pub.sh -M `pwd` $x >/dev/null   # to be improved
 fi
done

cd $outdir
framebatch_gapfill.sh -l -P -o 5 120 $rgl $azl
