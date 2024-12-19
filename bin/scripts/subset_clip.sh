#!/bin/bash

# this script can also work in full frame but it is set to allow subset clipping to even smaller area
# e.g. to find corner reflector where we want to load intensities into datacube

cd /gws/nopw/j04/nceo_geohazards_vol1/projects/LiCS/proc/current/subsets/kladno/095D
centerlon=14.00393208465662
centerlat=50.12668894637102
halflenrg=14
halflenazi=3
outdir=subset.kz4

mkdir -p $outdir/RSLC
m=`ls SLC | head -n 1`
h=`ls GEOC.meta*/*geo.hgt.tif | head -n 1`
hei=`python3 -c "import rioxarray;a=rioxarray.open_rasterio('"$h"'); aa=a.sel(x="$centerlon",y="$centerlat", method='nearest'); print(aa.values[0])"`

slc=SLC/$m/$m.slc
slcpar=SLC/$m/$m.slc.par
coord_to_sarpix $slcpar - - $centerlat $centerlon $hei | grep "SLC/MLI range, azimuth pixel (int)" > tmp.subset.centre
azipx=`cat tmp.subset.centre | rev | gawk {'print $1'} | rev`
rgpx=`cat tmp.subset.centre | rev | gawk {'print $2'} | rev`

let azi2=azipx+$halflenazi
let azi1=azipx-$halflenazi
let azidiff=azi2-azi1+1

let rg2=rgpx+$halflenrg
let rg1=rgpx-$halflenrg
let rgdiff=rg2-rg1+1

# TODO from here
echo "clipping reference epoch"
mkdir $outdir/RSLC/$m
SLC_copy $slc $slcpar $outdir/RSLC/$m/$m.rslc $outdir/RSLC/$m/$m.rslc.par - - $rg1 $rgdiff $azi1 $azidiff - - >/dev/null 2>/dev/null
multi_look $outdir/RSLC/$x/$x.rslc $outdir/RSLC/$x/$x.rslc.par $outdir/RSLC/$x/$x.rslc.mli $outdir/RSLC/$x/$x.rslc.mli.par $rgl $azl >/dev/null 2>/dev/null


if [ $process_rslcs == 1 ]; then
	echo "performing full clipping"
	# ok, now clip the mosaics

	for x in `ls RSLC | grep 20`; do 
	 if [ -f RSLC/$x/$x.rslc ]; then
	 if [ ! -d $outdir/RSLC/$x ]; then
	   echo "clipping "$x
	   mkdir -p $outdir/RSLC/$x
	   SLC_copy RSLC/$x/$x.rslc RSLC/$x/$x.rslc.par $outdir/RSLC/$x/$x.rslc $outdir/RSLC/$x/$x.rslc.par - - $rg1 $rgdiff $azi1 $azidiff - - >/dev/null 2>/dev/null
	   # no need for multilooking here?... 
	   #multi_look $outdir/RSLC/$x/$x.rslc $outdir/RSLC/$x/$x.rslc.par $outdir/RSLC/$x/$x.rslc.mli $outdir/RSLC/$x/$x.rslc.mli.par $rgl $azl >/dev/null 2>/dev/null
	   # create_geoctiffs_to_pub.sh -M `pwd` $x >/dev/null   # to be improved
	 fi
	 fi
	done
fi
