#!/usr/bin/env python

# Jan 2021 - Milan Lazecky

import os, glob
import numpy as np
from scipy import interpolate

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

tif = '20210107_20210119.geo.diff_pha.gacos.notides.landmasked.tif'
newRasterfn = 'temp_nnmasked.tif'
dataset = gdal.Open(tif, gdal.GA_ReadOnly)
band = dataset.GetRasterBand(1)
ifg = band.ReadAsArray()

#ifg[ifg==0] = np.nan
#array = np.ma.masked_invalid(ifg)

array = np.ma.masked_equal(ifg,0)
#to fill only 0, keep nans..

GD1[GD1==np.nan] = 0
originX, pixelWidth, b, originY, d, pixelHeight = dataset.GetGeoTransform() 
driver = gdal.GetDriverByName('GTiff')
cols = GD1.shape[1]
rows = GD1.shape[0]
GDT_dtype = gdal.GDT_Float32
outRaster = driver.Create(newRasterfn, cols, rows, 1, GDT_dtype)
outRaster.SetGeoTransform((originX, pixelWidth, 0, originY, 0, pixelHeight))
outband = outRaster.GetRasterBand(1)
outband.WriteArray(GD1)
prj=dataset.GetProjection()
outRasterSRS = osr.SpatialReference(wkt=prj)
outRaster.SetProjection(outRasterSRS.ExportToWkt())
outband.FlushCache()



def add_noise(ifg, coh = None, factor=10):
    '''
    Simulate (speckle-like) Gaussian noise
     - the higher the factor, the higher noise (try between 1-100 - factor 10 would have random_mean=0 and std=1)
     - the factor would also change by coherence (should be value 0-1), if provided (as np.array)
    '''
    #use e.g.:
    # ifg = np.fromfile('filtdiff.5',dtype=np.complex64)
    # coh = np.fromfile('filtcoh.5',dtype=np.float32)
    # # coh = coh/255
    # coh = coh.reshape(len,wid) ...
    if type(coh) != type(None):
        noise = factor * (1/coh) * np.random.standard_normal(size=ifg.shape)/10
    else:
        noise = factor * np.random.standard_normal(size=ifg.shape)/10
    outifg = ifg + ifg*noise
    return outifg


def 
import matplotlib.pyplot as plt

factor=2
ifg=pha.z.values
fig, (ax1, ax2, ax3) = plt.subplots(1, 3)
fig.suptitle('ifg comparison')
#ax3.imshow(np.angle(add_noise(ifg, coh = coh, factor=factor)), cmap='RdYlBu')
ax3.imshow(np.angle(add_noise(ifg, coh = coh.z.values, factor=factor)), cmap='RdYlBu')
ax2.imshow(np.angle(add_noise(ifg, coh = None, factor=factor)), cmap='RdYlBu')
ax1.imshow(np.angle(ifg), cmap='RdYlBu')
fig.show()


#use ncs:
import xarray as xr
pha = xr.open_dataset('temp.ifg.masked.nc')
coh = xr.open_dataset('temp.coh.nc')
outifg = add_noise(pha.z.values, coh.z.values, 10)

out = pha.copy(deep=True)
out.z.values = outifg
out.to_netcdf('temp.ifg.noised.nc')
#musim to pak wrapovat jeste
# gmt grdmath temp.ifg.noised.nc WRAP = temp.ifg.masked.nc 

outifg2 = add_noise(out.z.values, coh.z.values, 2)
out2 = pha.copy(deep=True)
out2.z.values = outifg2
out2.to_netcdf('temp.ifg.noised2.nc')
