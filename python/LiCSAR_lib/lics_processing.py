#!/usr/bin/env python3

################################################################################
# LiCSAR InSAR Processing functions
# by Milan Lazecky, 2021-2022, University of Leeds
#
# version: 1.0.0 (2022-11-03)
#
################################################################################
#Imports
################################################################################
import os, glob
import shutil
import subprocess

import xarray as xr
xr.set_options(keep_attrs=True)
import rioxarray

import numpy as np
from scipy import interpolate


def mm2rad_s1(inmm, rad2mm=False):
    """Converts from mm to radians (for Sentinel-1)
    """
    #speed_of_light = 299792458 #m/s
    speed_of_light = 299702547  #m/s ... in case of all-in-air
    radar_freq = 5.405e9  #for S1
    wavelength = speed_of_light/radar_freq #meter
    coef_r2m = wavelength/4/np.pi*1000 #rad -> mm,
    if rad2mm:
        # apologies for the inmm/outrad naming
        outrad = inmm*coef_r2m
    else:
        outrad = inmm/coef_r2m
    return outrad


def rad2mm_s1(inrad):
    return mm2rad_s1(inrad, rad2mm=True)


def ztd2sltd(icamsh5, ugeo, outif = None):
    ''' Function to convert output of ICAMS ZTD to SLTD (in corresponding coordinates)'''
    ugeo = rioxarray.open_rasterio(ugeo)
    icams = xr.open_dataset(icamsh5)
    icams = icams.interp_like(ugeo, method='linear')
    icams = icams.where(~np.nan(ugeo))
    icams = icams.where(ugeo != 0) # just in case..
    icams = icams/ugeo
    if outif:
        ugeo.values = icams.values
        ugeo.rio.to_raster(outif, compress='deflate')
    return icams


def interpolate_nans(array, method='cubic'):
    """Interpolation of NaN values in a grid

    Args:
        array (np.array): numpy array with nans to interpolate
        method (string): interpolation method for griddata function, e.g. cubic

    Returns:
        np.array: interpolated grid
    """
    array = np.ma.masked_invalid(array)
    x = np.arange(0, array.shape[1])
    y = np.arange(0, array.shape[0])
    xx, yy = np.meshgrid(x, y)
    x1 = xx[~array.mask]
    y1 = yy[~array.mask]
    newarr = array[~array.mask]
    GD1 = interpolate.griddata((x1, y1), newarr.ravel(),(xx, yy),method=method)
    GD1 = np.array(GD1)
    return GD1


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


def load_tif2xr(tif, cliparea_geo=None, tolonlat=True):
    """loads geotiff to xarray.DataArray
    
    Args:
        tif (string): path to geotiff
        cliparea_geo (string): use GMT/LiCSBAS string to identify area to clip, in geo-coordinates, as ``'lon1/lon2/lat1/lat2'``
        tolonlat (boolean): if True, return as lon lat coordinates
    
    Returns:
        xr.DataArray: loaded contents
    """
    xrpha = rioxarray.open_rasterio(tif)
    xrpha = xrpha.squeeze('band')
    xrpha = xrpha.drop('band')
    
    if cliparea_geo:
        minclipx, maxclipx, minclipy, maxclipy = cliparea_geo2coords(cliparea_geo)
        xrpha = xrpha.sel(x=slice(minclipx, maxclipx), y=slice(maxclipy, minclipy))
    if tolonlat:
        xrpha = xrpha.rename({'x': 'lon','y': 'lat'})
    return xrpha


def cliparea_geo2coords(cliparea_geo):
    """Exports the string to min/max clip values

    Args:
        cliparea_geo (str): clip boundaries, e.g. 'lon1/lon2/lat1/lat2'

    Returns:
        float, float, float, float: minclipx, maxclipx, minclipy, maxclipy
    """
    minclipx, maxclipx, minclipy, maxclipy = cliparea_geo.split('/')
    minclipx, maxclipx, minclipy, maxclipy = float(minclipx), float(maxclipx), float(minclipy), float(maxclipy)
    if minclipy > maxclipy:
        print('you switched min max in crop coordinates (latitude). fixing')
        tmpcl = minclipy
        minclipy = maxclipy
        maxclipy = tmpcl
    if minclipx > maxclipx:
        print('you switched min max in crop coordinates (longitude). fixing')
        tmpcl = minclipx
        minclipx = maxclipx
        maxclipx = tmpcl
    return minclipx, maxclipx, minclipy, maxclipy


