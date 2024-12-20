#!/bin/bash

# this script can also work in full frame but it is set to allow subset clipping to even smaller area
# e.g. to find corner reflector where we want to load intensities into datacube

makemlinc=1  # this will create MLI NetCDF file.. that should be actually converted to amplitude (and put dB in)

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

if [ $makemlinc == 1 ]; then
# 1. generate the small size mlis (lazy approach. not effective. still ok)
for r in `ls RSLC`; do
 # TODO from here
 #echo "clipping reference epoch"
 echo $r
 mkdir -p $outdir/RSLC/$r
 if [ ! -f $outdir/RSLC/$r/$r.rslc ]; then
  SLC_copy RSLC/$r/$r.rslc RSLC/$r/$r.rslc.par $outdir/RSLC/$r/$r.rslc $outdir/RSLC/$r/$r.rslc.par - - $rg1 $rgdiff $azi1 $azidiff - - >/dev/null 2>/dev/null
 fi
 if [ ! -f $outdir/RSLC/$r/$r.rslc.mli ]; then
  multi_look $outdir/RSLC/$r/$r.rslc $outdir/RSLC/$r/$r.rslc.par $outdir/RSLC/$r/$r.rslc.mli $outdir/RSLC/$r/$r.rslc.mli.par 1 1 >/dev/null 2>/dev/null
 fi
done
# now load that to netcdf (and convert to amplitude)

cd $outdir
# python...
import os
import numpy as np
import xarray as xr
import datetime as dt

width=$rgdiff # 29
length=$azidiff # 7
# for track 95:
time95 = '05 18 2.04720'
h = 5
m = 18
s = int(2.042720)

epochs = []
data = []
for epoch in os.listdir('RSLC'):
    print(epoch)
    mli='RSLC/'+epoch+'/'+epoch+'.rslc.mli'
    if os.path.exists(mli):
        a=np.fromfile(mli, dtype=np.float32)
        ampdb=10*np.log10(np.sqrt(a.byteswap())).reshape((length, width))
        epochdt = dt.datetime(int(epoch[:4]), int(epoch[4:6]), int(epoch[6:]), h, m, s)
        epochs.append(epochdt)
        data.append(ampdb)

bb=xr.DataArray(data=data, dims=['epoch','azi','rg'], coords=dict(azi=np.arange(length), rg=np.arange(width), epoch=epochs))
bb=xr.Dataset(data_vars=dict(amplitude_dB=bb))
outnc = '/gws/nopw/j04/nceo_geohazards_vol1/public/shared/temp/earmla/kladno/kz4/095.nc'
bb.to_netcdf(outnc)



fi
