#!/usr/bin/env python3

################################################################################
# Export of LiCSAR data towRDS STAMPS
# by Milan Lazecky, 2021-2022, University of Leeds
#
# version: 1.0.0 (2022-06-24)

# STARTING PARAMETERS:

import os
import pandas as pd
from framecare import get_master
import xarray as xr
import rioxarray
import numpy as np
import subprocess as subp

frame='062D_06031_131313'
track=str(int(frame[:3]))
geoframedir=os.path.join(os.environ['LiCSAR_public'],track,frame)
cliparea='102.5/103.4/30.1/30.7'


#outputs

outdir='stampsproc'
if not os.path.exists(outdir):
    os.mkdir(outdir)

outslcdir=os.path.join(outdir,'slc')
if not os.path.exists(outslcdir):
    os.mkdir(outslcdir)

outgeodir=os.path.join(outdir,'geo')
if not os.path.exists(outgeodir):
    os.mkdir(outgeodir)

outifgdir=os.path.join(outdir,'diff0')
if not os.path.exists(outifgdir):
    os.mkdir(outifgdir)





def extract_lonlat(xrda, outdir):
    lonfile=os.path.join(outdir,'lon')
    latfile=os.path.join(outdir,'lat')
    if os.path.exists(lonfile):
        return False
    dtype=np.float32
    #xrda.lon.values.astype(dtype).tofile(lonfile)
    #xrda.lat.values.astype(dtype).tofile(latfile)
    np.tile(xrda.lon.values.astype(dtype), (len(xrda.lat), 1)).tofile(lonfile) #.shape
    np.tile(xrda.lat.values.astype(dtype), (len(xrda.lon), 1)).T.tofile(latfile) #.shape
    return True


def get_bperp(pair,baselines):
    bperp=0
    for epoch in pair.split('_'):
        if int(epoch) in baselines.index:
            bperp1=baselines.loc[int(epoch)].bperp
        else:
            print('warning - epoch {} not in baselines file. using zero'.format(epoch))
            bperp1=0
        bperp=bperp1-bperp
    return bperp


# this can be taken as:
# from lics_unwrap import load_tif2xr, magpha2RI_array

def load_tif2xr(tif, cliparea_geo=None, tolonlat=True):
    """loads geotiff to xarray.DataArray
    
    Args:
        tif (string): path to geotiff
        cliparea_geo (string): use GMT/LiCSBAS string to identify area to clip, in geo-coordinates, as 'lon1/lon2/lat1/lat2'
        tolonlat (boolean): if True, return as lon lat coordinates
    
    Returns:
        xr.DataArray: loaded contents
    """
    xrpha = rioxarray.open_rasterio(tif)
    xrpha = xrpha.squeeze('band')
    xrpha = xrpha.drop('band')
    
    if cliparea_geo:
        minclipx, maxclipx, minclipy, maxclipy = cliparea_geo.split('/')
        minclipx, maxclipx, minclipy, maxclipy = float(minclipx), float(maxclipx), float(minclipy), float(maxclipy)
        xrpha = xrpha.sel(x=slice(minclipx, maxclipx), y=slice(maxclipy, minclipy))
    if tolonlat:
        xrpha = xrpha.rename({'x': 'lon','y': 'lat'})
    return xrpha


def magpha2RI_array(mag, pha):
    """Converts arrays of magnitude and phase to complex number array (real and imaginary)

    Args:
        mag (np.array): numpy array with magnitude values
        pha (np.array): numpy array with phase values

    Returns:
        np.array: complex number array
    """
    R = np.cos(pha) * mag
    I = np.sin(pha) * mag
    out = R + 1j*I
    return out


def get_param_par(mlipar, field):
    """
    Get parameter from mli.par or dem_par file. Examples of fields are;
     - range_samples
     - azimuth_lines
     - range_looks
     - azimuth_looks
     - range_pixel_spacing (m)
     - azimuth_pixel_spacing (m)
     - radar_frequency  (Hz)
    """
    value = subp.check_output(['grep', field,mlipar]).decode().split()[1].strip()
    return value


def get_spacing_xr(ifg):
    '''Get spacing of given xarray dataset (must have lon, lat coords). Assuming rectangle spacing
    '''
    resdeg = np.abs(ifg.lat[1]-ifg.lat[0])
    latres = 111.32 * np.cos(np.radians(ifg.lat.mean()))
    res = float(latres * resdeg)*1000 # in m
    return res


# metadata and hgt:
hgtfile = os.path.join(geoframedir,'metadata', frame+'.geo.hgt.tif')
outhgt =  os.path.join(outgeodir,'dem')
hgt = load_tif2xr(hgtfile, cliparea_geo=cliparea)
hgt.values.tofile(outhgt)
rc=extract_lonlat(hgt, outgeodir)

length=len(hgt.lat)
width=len(hgt.lon)
cmd='echo '+str(len(hgt.lat))+' > '+os.path.join(outdir,'len.txt')
rc=os.system(cmd)
cmd='echo '+str(len(hgt.lon))+' > '+os.path.join(outdir,'width.txt')
rc=os.system(cmd)

basefile=os.path.join(geoframedir,'metadata','baselines')
cmd='cp '+basefile+' '+outdir
rc=os.system(cmd)

metafile=os.path.join(geoframedir,'metadata','metadata.txt')
cmd='cp '+metafile+' '+outdir
rc=os.system(cmd)







print('Extracting information from/to reference mli par file')
rdframedir=os.path.join(os.environ['LiCSAR_procdir'],track,frame)
slcs=os.listdir(os.path.join(rdframedir,'SLC'))
master=get_master(frame)
mlipar=os.path.join(rdframedir,'SLC',master,master+'.slc.mli.par')
parfile=os.path.join(outslcdir,os.path.basename(mlipar))

fieldsfrompar = ['image_format',
'near_range_slc',
'sar_to_earth_center',
'earth_radius_below_sensor',
'center_range_slc',
'prf',
'heading',
'sensor']

f = open(parfile, 'w')
for field in fieldsfrompar:
    value = get_param_par(mlipar, field)
    f.write(field+':     '+value+'\n')

# wid and len:
field='azimuth_lines'
value=str(length)
f.write(field+':     '+value+'\n')
field='range_samples'
value=str(width)
f.write(field+':     '+value+'\n')

# resolution in 'm':
spacing = get_spacing_xr(hgt)
field='range_pixel_spacing'
ovalue_rgs = get_param_par(mlipar, field)
value = str(spacing)
f.write(field+':     '+value+'\n')
field='azimuth_pixel_spacing'
ovalue_azs = get_param_par(mlipar, field)
f.write(field+':     '+value+'\n')

# looks:
field='range_looks'
ovalue_rgl = get_param_par(mlipar, field)
orig_spacing_rg = float(ovalue_rgs)/float(ovalue_rgl)
value = round(spacing/orig_spacing_rg)
f.write(field+':     '+str(value)+'\n')
field='azimuth_looks'
ovalue_azl = get_param_par(mlipar, field)
orig_spacing_az = float(ovalue_azs)/float(ovalue_azl)
value = round(spacing/orig_spacing_az)
f.write(field+':     '+str(value)+'\n')

f.close()



'''
# we decided we don't need MLIs
print('extracting reference MLI')
mlitif= os.path.join(geoframedir,'epochs',master,master+'.geo.mli.tif')
outmli=os.path.join(outslcdir,master+'.mli')
if os.path.exists(mlitif) and not os.path.exists(outmli):
    mli = load_tif2xr(mlitif, cliparea_geo=cliparea)
    if mli.shape != hgt.shape:
        print('wrong size of MLI. skipping')
    else:
        mli.astype(np.float32).values.tofile(outmli)
'''

pairset = os.listdir(os.path.join(geoframedir,'interferograms')) 
print('extracting ifgs (warning, should use glob instead of listdir, might cause errors)')
# also prepare baselines:
baselines = pd.read_csv(basefile, header=None, delimiter=' ')
baselines = baselines[[1,2]]
baselines = baselines.drop_duplicates()
baselines=baselines.rename(columns={1: "epoch", 2: "bperp"})
baselines=baselines.set_index('epoch')

for pair in pairset:
    #pair='20220524_20220605'
    geoifgdir = os.path.join(geoframedir,'interferograms',pair)
    difftif = os.path.join(geoifgdir,pair+'.geo.diff_unfiltered_pha.tif')
    cohtif = os.path.join(geoifgdir,pair+'.geo.cc.tif')
    outifg = os.path.join(outifgdir,pair,pair+'.diff')
    outcc = os.path.join(outifgdir,pair,pair+'.cc')
    outbperp = os.path.join(outifgdir,pair,pair+'.bperp')
    if os.path.exists(difftif) and os.path.exists(cohtif) and not os.path.exists(outifg):
        if not os.path.exists(os.path.join(outifgdir,pair)):
            os.mkdir(os.path.join(outifgdir,pair))
        pha = load_tif2xr(difftif, cliparea_geo=cliparea)
        if pha.shape != hgt.shape:
            print('wrong size of ifg '+pair+'. skipping')
            continue
        #print('extracting ifg '+pair)
        coh = load_tif2xr(cohtif, cliparea_geo=cliparea)
        coh = coh/255
        cpx = magpha2RI_array(coh, pha)
        cpx.astype(np.complex64).values.tofile(outifg)
        coh.astype(np.float32).values.tofile(outcc)
        bperp=get_bperp(pair,baselines)
        with open(outbperp, 'w') as the_file:
            the_file.write(str(bperp))