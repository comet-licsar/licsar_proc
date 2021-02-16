#!/bin/bash
source $LiCSARpath/lib/LiCSAR_bash_lib.sh
#10 is approx coh 0.04
cohthr=10



#multiscale here:
inpha=tin1.nc  #this is masked nc
maskin=mask.fullin.nc
fullmask=mask.nc
cohin=coh.nc #0-1
for scale in 5 3 2; do
 
#downsample both pha and coh, and mask... in 5x5
cat << EOF > snaphu.conf
STATCOSTMODE  DEFO
INFILEFORMAT  FLOAT_DATA
CORRFILEFORMAT  FLOAT_DATA
OUTFILEFORMAT FLOAT_DATA
ESTFILEFORMAT FLOAT_DATA
RMTMPTILE TRUE
EOF


#time snaphu -f snaphu.conf -M mask1 -o unw1 -c coh1 -g unw1.conncomp pha1 $width #-e est1

#for snaphu, the bytemask is in signed char
import xarray as xr
import numpy as np
import os
import cv2
from scipy import interpolate

inphaf='tin1.nc'
maskinf='mask.fullin.nc'
cohinf='coh.nc'
fullmaskf='mask.nc'

def interpolate_nans(array):
    array = np.ma.masked_invalid(array)
    x = np.arange(0, array.shape[1])
    y = np.arange(0, array.shape[0])
    xx, yy = np.meshgrid(x, y)
    x1 = xx[~array.mask]
    y1 = yy[~array.mask]
    newarr = array[~array.mask]
    GD1 = interpolate.griddata((x1, y1), newarr.ravel(),(xx, yy),method='nearest') #, fill_value=0)
    GD1 = np.array(GD1)
    return GD1


inpha = xr.open_dataset(inphaf)
incoh = xr.open_dataset(cohinf)
inmask = xr.open_dataset(maskinf)
fullmask = xr.open_dataset(fullmaskf)

prevscale = None
for scale in [5,3,2,1]:
    pha_coarsened = inpha.coarsen({'lat': scale, 'lon': scale}, boundary='pad').median()
    coh_coarsened = incoh.coarsen({'lat': scale, 'lon': scale}, boundary='pad').median()
    mask_coarsened = inmask.coarsen({'lat': scale, 'lon': scale}, boundary='pad').median()
    fullmask_coarsened = fullmask.coarsen({'lat': scale, 'lon': scale}, boundary='pad').median()
    #wrap it and shift by pi, as snaphu expects:
    pha_coarsened = pha_coarsened.fillna(10)
    pha_masked = pha_coarsened.copy(deep=True)
    pha_masked.z.values = pha_masked.z.values * mask_coarsened.z.values
    pha_masked = pha_masked.where(pha_masked.z != 0)
    #now use griddata, finally!
    pha_masked.z.values = interpolate_nans(pha_masked.z.values)
    ######pha_masked.to_netcdf('pha'+str(scale)+'.nc')
    #now replace those 10s to zeroes..
    pha_masked = pha_masked.where(pha_masked.z != 10)#.fillna(0)
    phases = np.array(pha_masked.z.values).astype(np.float)
    pha_masked.z.values = np.arctan2(np.sin(phases), np.cos(phases)) + np.pi
    pha_masked.z.values = pha_masked.z.fillna(0).values
    #phases.tofile('pha'+str(scale))
    pha_masked.to_netcdf('pha'+str(scale)+'.nc')
    coh_coarsened.to_netcdf('coh'+str(scale)+'.nc')
    mask_coarsened.to_netcdf('mask'+str(scale)+'.nc')
    rc = os.system('gmt grd2xyz -ZTLf -bof pha{0}.nc > pha{0}'.format(str(scale)))
    rc = os.system('gmt grd2xyz -ZTLf -bof coh{0}.nc > coh{0}'.format(str(scale)))
    rc = os.system('gmt grd2xyz -ZTLc -bof mask{0}.nc > mask{0}'.format(str(scale)))
    #os.system('gmt grdfill temppha5.nc -An -Gtemppha5.filtered.nc')
    #put nans instead of 0 (=masked)
    #pha_masked.interpolate_na(method='nearest')
    #finally export to binaries
    #np.array(pha_masked.z).astype(np.float).tofile('pha'+str(scale))
    #np.array(coh_coarsened.z).astype(np.float).tofile('coh'+str(scale))
    #np.array(mask_coarsened.z).astype(np.byte).tofile('mask'+str(scale))
    
    
    #and now start snaphu processing
    width = len(pha_coarsened.lon)
    
    if prevscale:
        # i should load the prevscale unw and interpolate it, as est$scale !
        dimopposite = (pha_coarsened.z.values.shape[1], pha_coarsened.z.values.shape[0] )
        est = cv2.resize(unwdone,dsize=dimopposite, interpolation=cv2.INTER_CUBIC)
        estxr = pha_coarsened.copy(deep=True)
        estxr.z.values = est - prevmedian #np.nanmedian(estxr.z.values)
        estxr.to_netcdf('est{}.nc'.format(scale))
        rc = os.system('gmt grd2xyz -ZTLf -bof est{0}.nc > est{0}'.format(str(scale)))
        #tempmask = xr.open_dataset('mask{}.nc'.format(prevscale))
        rc = os.system('time snaphu -f snaphu.conf -M mask{0} -o unw{0} -c coh{0} -e est{0} pha{0} {1}'.format(str(scale),str(width)))
    else:
        #first scaling = no estimates
        rc = os.system('time snaphu -f snaphu.conf -M mask{0} -o unw{0} -c coh{0} pha{0} {1} #-e est1'.format(str(scale),str(width)))
    
    # and now import and save as nc
    unwdone = np.fromfile('unw{}'.format(scale),dtype=np.float32)
    a = fullmask_coarsened.copy(deep=True)
    unwdone = unwdone.reshape(fullmask_coarsened.z.shape)
    unwdone = np.flip(unwdone,axis=0)
#
    a.z.values = unwdone*fullmask_coarsened.z.values
    a.z.values[a.z.values==0] = np.nan
    prevmedian = np.nanmedian(a.z.values)
    #great, now make it ready for the next step
    #est = fullmask_coarsened.copy(deep=True)
    prevscale = scale
    prevunw = unwdone.copy()
    
    #and just save the unw to nc
    a.z.values = a.z.values - np.nanmedian(a.z.values)
    a.to_netcdf('unw{}.nc'.format(scale))






for scale in 5 3 2; do
 gmt grdmath -N pha$scale.nc 0 NAN 10 DENAN mask$scale.nc MUL 0 NAN = temp.tofill$scale.nc


    os.system('')
    phases.tofile('pha'+str(scale))
    np.array(coh_coarsened.z).astype(np.float).tofile('coh'+str(scale))
    np.array(mask_coarsened.z).astype(np.byte).tofile('mask'+str(scale))
    width = len(inpha.lon)
    os.system('time snaphu -f snaphu.conf.multiscale {}'.format(width))
# dont forget pha must be in 0-2pi!
#convert to snaphu binaries and process them by snaphu
time snaphu -f snaphu.conf.3.assisted $width
snaphu.conf.1.assisted.10to5casc

#then interpolate them to lower step
done



def upsize(a, wid):
est = cv2.resize(unwdone,dsize=(origwid,origlen), interpolation=cv2.INTER_CUBIC)
>>> import cv2                                                                                                            
>>> import numpy as np                                                                                                    
>>> import pandas as pd                                                                                                   
>>> import os                                                                                                             
>>> wid = pd.read_csv('wid.3',header=None)                                                                                
>>> wid = int(wid[0])                                                                                                                                                          
>>> unwest = np.fromfile('unw10to5to3.3',dtype=np.float32)                                                                
>>> length = int(len(unwest)/wid)                            
>>> unwest = unwest.reshape(length,wid)                                                                                   
>>> origwid = pd.read_csv('wid.txt',header=None)                                                                          
>>> origwid = int(origwid[0])                                                                                             
>>> fsize = os.path.getsize('filtdiff')                                                                                   
>>> origlen = int(fsize/8/origwid)                                                                                        
>>> a = cv2.resize(unwest,dsize=(origwid,origlen), interpolation=cv2.INTER_CUBIC)                                         
>>> a = a - np.nanmedian(a)                                                                                               
>>> a.tofile('unw10to5to3to1.cv2')                                                                                        
>>> exit()  



# M. Lazecky, 2021

if [ -z $2 ]; then
 echo "parameters: frame ifg"
 exit
fi
frame=$1
ifgid=$2

track=`track_from_frame $frame`
heredir=`pwd`
#this script will unwrap geocoded data...
#just set those files:
maskfile=$LiCSAR_public/$track/$frame/metadata/$frame.geo.landmask.tif
ifg=$LiCSAR_public/$track/$frame/interferograms/$ifgid/$ifgid.geo.diff_pha.tif
coh=$LiCSAR_public/$track/$frame/interferograms/$ifgid/$ifgid.geo.cc.tif
outunw=$LiCSAR_public/$track/$frame/interferograms/$ifgid/$ifgid.geo.unw.tif

cd $LiCSAR_public/$track/$frame/interferograms/$ifgid
mkdir temp 2>/dev/null

width=`gmt grdinfo $ifg | grep n_columns | rev | gawk {'print $1'} | rev`

#preparing mask.nc
echo "preparing masks"
gmt grdcut -N1 $maskfile -Gtemp/mask.landmask.nc -R$ifg  #grid reg
gmt grdedit temp/mask.landmask.nc -T -R$ifg   #to pixel reg
gmt grdmath $coh 0 NAN 0 GT 0 DENAN = temp/mask.outarea.nc
gmt grdmath $coh 0 NAN $cohthr GT 1 DENAN = temp/mask.inarea.nc #pixel reg!!!! 
#gmt grdedit temp/mask.inarea.nc -T -R$ifg   #to pixel reg
gmt grdmath -N temp/mask.inarea.nc temp/mask.landmask.nc MUL = temp/mask.fullin.nc
gmt grdmath -N temp/mask.outarea.nc temp/mask.fullin.nc MUL = temp/mask.nc

#gapfilling masked areas
echo "gapfilling masked areas"
#cp temp/mask.nc temp/mask.georeg.nc
#
# ok, masking only the internal parts..
cp temp/mask.fullin.nc temp/mask.fullin.georeg.nc
#gmt grdedit temp/mask.georeg.nc -T -R$ifg   #to pixel reg
gmt grdedit temp/mask.fullin.georeg.nc -T -R$ifg   #to pixel reg
#gmt grdmath -N $ifg temp/mask.fullin.georeg.nc MUL 0 NAN = temp/ifg.fullin.tofill.nc



#for ifgnc in tin1.nc tin2.nc tin3.nc tin5.nc tin10.nc; do
for ifgnc in tin2.nc tin3.nc tin5.nc; do
    gmt grdmath -N $ifgnc 0 NAN 10 DENAN temp/mask.fullin.georeg.nc MUL 0 NAN = temp/ifg.masked.tofill.nc

    #gmt grdmath -N $ifg temp/mask.georeg.nc MUL 0 NAN = temp/ifg.masked.nc
    #time gmt grdfill temp/ifg.masked.nc -An -Gtemp/pha1.nc
    gmt grdfill temp/ifg.masked.tofill.nc -An -Gtemp/pha1.filled.nc
    gmt grdmath temp/pha1.filled.nc 10 NAN 0 NAN 0 DENAN = temp/pha1.nc

    #make binaries for snaphu:
    echo "unwrapping by snaphu"
    #input is phase only = snaphu needs phase in 0,2 pi
    gmt grdmath temp/pha1.nc PI ADD = temp/pha1.02pi.nc
    gmt grd2xyz -ZTLf -bof temp/pha1.02pi.nc > temp/pha1
    
    cd temp
    time snaphu -f snaphu.conf -M mask1 -o unw1 -c coh1 -g unw1.conncomp pha1 $width #-e est1
    #time snaphu -f snaphu.conf -M mask1 -o unw1 -c coh1 pha1 $width #-e est1
#now convert back to nc and re-include mask
echo "converting to geotiff"

    python3 unw2nc.py 
    mv unw1.nc unw_$ifgnc
    #make preview
    create_preview_unwrapped unw_$ifgnc 
    cp unw1.conncomp $ifgnc.conncomp
    display unw_`echo $ifgnc | cut -d '.' -f1`png &
    #gmt grdconvert -G$outunw=gd:GTiff unw1.nc
    #mv unw1png `echo $outunw | rev | cut -c 4- | rev`png
    cd ..;# rm -r temp
done
#cd $heredir
#done
#echo "now take a look:"
#echo "display unw1.png"

