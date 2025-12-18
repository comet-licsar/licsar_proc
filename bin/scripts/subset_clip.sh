#!/bin/bash

# this script can also work in full frame but it is set to allow subset clipping to even smaller area
# e.g. to find corner reflector where we want to load intensities into datacube
# usage with params:
# subset_clip.sh lon lat locationname [half range] [half azi]
# outputs will be stored as subset.$locationname
#
# must be run inside the subset directory, e.g. in
# /gws/nopw/j04/nceo_geohazards_vol1/projects/LiCS/proc/current/subsets/kladno/095D
makemlinc=1  # this will create MLI NetCDF file.. that should be actually converted to amplitude (and put dB in)

minepoch=20240101

#ls /gws/nopw/j04/nceo_geohazards_vol1/public/shared/temp/earmla/kladno/kz4
#cd /gws/nopw/j04/nceo_geohazards_vol1/projects/LiCS/proc/current/subsets/kladno/095D
centerlon=14.00393208465662
centerlat=50.12668894637102
halflenrg=14
halflenazi=3
locname=kz4

# let's get it from params:
if [ ! -z $3 ]; then
centerlon=$1
centerlat=$2
locname=$3
if [ ! -z $5 ]; then
  halflenrg=$4
  halflenazi=$5
fi
else
  echo " please provide parameters"
  exit
fi
outdir=subset.$locname
relorb=`pwd | rev | cut -d '/' -f 1 | rev`

mkdir -p $outdir/RSLC
cp $outdir/command.in $outdir/command.in.old 2>/dev/null
echo "cd "`pwd`"; subset_clip.sh "$centerlon" "$centerlat" "$locname" "$halflenrg" "$halflenazi > $outdir/command.in
m=`ls SLC | head -n 1`
h=`ls GEOC.meta*/*geo.hgt.tif | head -n 1`

slc=SLC/$m/$m.slc
slcpar=SLC/$m/$m.slc.par
if [ ! -f $outdir/tmp.subset.centre ]; then
 #hei=`python3 -c "import rioxarray;a=rioxarray.open_rasterio('"$h"'); aa=a.sel(x="$centerlon",y="$centerlat", method='nearest'); print(aa.values[0])"`
 hei=`gdallocationinfo -geoloc -valonly $h $centerlon $centerlat`
 coord_to_sarpix $slcpar - - $centerlat $centerlon $hei | grep "SLC/MLI range, azimuth pixel (int)" > $outdir/tmp.subset.centre
fi
azipx=`cat $outdir/tmp.subset.centre | rev | gawk {'print $1'} | rev`
rgpx=`cat $outdir/tmp.subset.centre | rev | gawk {'print $2'} | rev`

let azi2=azipx+$halflenazi
let azi1=azipx-$halflenazi
let azidiff=azi2-azi1+1

let rg2=rgpx+$halflenrg
let rg1=rgpx-$halflenrg
let rgdiff=rg2-rg1+1


# 1. generate the small size mlis (lazy approach. not effective. still ok)
echo "clipping rslcs"
for r in `ls RSLC`; do
 if [ $r -ge $minepoch ]; then
 # TODO from here
 #echo "clipping reference epoch"
 echo $r
 mkdir -p $outdir/RSLC/$r
 if [ ! -f $outdir/RSLC/$r/$r.rslc ]; then
  SLC_copy RSLC/$r/$r.rslc RSLC/$r/$r.rslc.par $outdir/RSLC/$r/$r.rslc $outdir/RSLC/$r/$r.rslc.par - - $rg1 $rgdiff $azi1 $azidiff - - >/dev/null 2>/dev/null
 fi
 fi
done

if [ $makemlinc == 1 ]; then
  echo "generating mlis"
for r in `ls RSLC`; do
   if [ $r -ge $minepoch ]; then
    echo $r
 if [ ! -f $outdir/RSLC/$r/$r.rslc.mli ]; then
  multi_look $outdir/RSLC/$r/$r.rslc $outdir/RSLC/$r/$r.rslc.par $outdir/RSLC/$r/$r.rslc.mli $outdir/RSLC/$r/$r.rslc.mli.par 1 1 >/dev/null 2>/dev/null
 fi
   fi
done
# radcal
  echo "basic radiometric calibration of mlis"
for r in `ls RSLC`; do
   if [ $r -ge $minepoch ]; then
     echo $r
 if [ ! -f $outdir/RSLC/$r/$r.rslc.mli.calibrated ]; then
  radcal_MLI $outdir/RSLC/$r/$r.rslc.mli $outdir/RSLC/$r/$r.rslc.mli.par - $outdir/RSLC/$r/$r.rslc.mli.calibrated - 1 >/dev/null 2>/dev/null
 fi
   fi
done

# now load that to netcdf (and convert to amplitude)
cat << EOF >$outdir/tonc.py
import os
import numpy as np
import xarray as xr
import datetime as dt
import LiCSAR_misc as misc

minepoch='20240101'

width=int(1+2*$halflenrg)  # 29 #$rgdiff # 29
length=int(1+2*$halflenazi)  #7 #$azidiff # 7
# for track 95:
# nano subset.kz4/RSLC/$m/$m.rslc.mli.par

def str2dt(strtime):
    #strtime = ['2024', '09', '08', '05', '09', '48.65361']
    # Convert to appropriate types
    year, month, day = map(int, strtime[:3])
    hour, minute = map(int, strtime[3:5])
    #second = float(strtime[5])
    second = int(strtime[5].split('.')[0])
    # Create datetime object
    return dt.datetime(year, month, day, hour, minute, second) #int(second), int((second % 1) * 1e6))


def mlipartime(mlipar):
    strtime=misc.grep1line('date',mlipar).split()[1:]
    return str2dt(strtime)


epochs = []
data = []
todb = False
for epoch in os.listdir('RSLC'):
    if (int(epoch)>int(minepoch)):
        print(epoch)
        mli='RSLC/'+epoch+'/'+epoch+'.rslc.mli.calibrated'
        mlipar='RSLC/'+epoch+'/'+epoch+'.rslc.mli.par'
        if not os.path.exists(mli):
            mli='RSLC/'+epoch+'/'+epoch+'.rslc.mli'
        if os.path.exists(mli):
            print(mli)
            a=np.fromfile(mli, dtype=np.float32)
            if todb:
                amp=10*np.log10(np.sqrt(a.byteswap())).reshape((length, width))
            else:
                amp=np.sqrt(a.byteswap()).reshape((length, width))
            # amp=np.flipud(amp) # should be in place..
            #if not np.isinf(ampdb):
            epochdt = mlipartime(mlipar)
            # epochdt = dt.datetime(int(epoch[:4]), int(epoch[4:6]), int(epoch[6:]), h, m, s)
            epochs.append(epochdt)
            data.append(amp)


bb=xr.DataArray(data=data, dims=['epoch','azi','rg'], coords=dict(azi=np.arange(length), rg=np.arange(width), epoch=epochs))
bb=xr.Dataset(data_vars=dict(amplitude=bb))
outnc = '$locname.$relorb.amp.nc'
# bb=xr.Dataset(data_vars=dict(amplitude_dB=bb))
#outnc = '/gws/nopw/j04/nceo_geohazards_vol1/public/shared/temp/earmla/kladno/kz4/095c.nc'
#outnc = '/gws/nopw/j04/nceo_geohazards_vol1/public/shared/temp/earmla/kladno/kz4/146.nc'
#outnc = '/gws/nopw/j04/nceo_geohazards_vol1/public/shared/temp/earmla/kladno/kz4/44.nc'
bb.to_netcdf(outnc)

EOF
echo "converting to NetCDF"
cd $outdir; python3 tonc.py

fi


exit


# how i did kladno:
cd /gws/nopw/j04/nceo_geohazards_vol1/projects/LiCS/proc/current/subsets/kladno
python3
import pandas as pd
import os, glob
csv='odrazece.csv'
df=pd.read_csv(csv, index_col=0)
for relorb in glob.glob('*[A,D]'):
    relorbpath = os.path.join(os.getcwd(), relorb)
    for i,row in df.iterrows():
        cmd="cd {0}; subset_clip.sh {1} {2} {3}".format(relorbpath, str(row['lon']), str(row['lat']), row['misto'])
        print(cmd)
        os.system(cmd)


# oh and... how to convert this to RCS [dBm^2] for direct comparison:
#For a point target like a corner reflector, the radar equation gives:
#\sigma = \sigma^0 \cdot A_{\text{res}}
#Where:
#- \sigma is the Radar Cross Section in m²
#- \sigma^0 is the normalized backscatter from the calibrated image
#- A_{\text{res}} is the resolution cell area in m²:
#A_{\text{res}} = \Delta r \cdot \Delta a
#- \Delta r: range resolution (m)
#- \Delta a: azimuth resolution (m)

import xarray as xr
import numpy as np
import glob,os
rngres = 2.3
azires = 14
for nc in glob.glob('*/subset.*/*.amp.nc'):
    print(nc)
    # nc='146A/subset.KZ-4/KZ-4.146A.amp.nc'
    outnc = nc.replace('.amp.nc','.rcs.nc')
    if os.path.exists(outnc):
        os.remove(outnc)
    nc=xr.open_dataset(nc)
    # \sigma_{\text{dBm}^2} = 10 \cdot \log_{10}(\sigma) + 30
    nc['RCS_dB'] = nc['amplitude'].copy()
    nc['RCS_dB'].values = 10*np.log10(nc.amplitude.where(nc.amplitude!=0).values*rngres*azires) #+30 # ok, unit will be dB, not dBm^2
    nc[['RCS_dB']].to_netcdf(outnc)


# finally, to identify the reflector based on max value in the last 5 images, and get the x,y to plot this..
import xarray as xr
import numpy as np
import glob
ncs = glob.glob('Tu-2.*')
nc = ncs[0]  # jako ukazka...
nc = xr.open_dataset(nc)
# prumerny RCS z poslednich 5 snimku:
avgrcs = nc['RCS_dB'][-5:].mean(axis=0)
# najit koordinaty max hodnoty:
maxazi, maxrg = np.unravel_index(avgrcs.argmax().item(), avgrcs.shape)
toplot = nc.sel(rg=maxrg, azi=maxazi)['RCS_dB']
# ted muzes treba:
toplot.plot()
# anebo vytahnout
x = toplot.epoch.values
y = toplot.values
