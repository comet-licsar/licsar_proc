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
from affine import Affine
import numpy as np
from scipy import interpolate
import statsmodels.api as sm

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
    # converting icams to standard geofile -- expecting tot_delay variable..
    sizey, sizex = icams.tot_delay.shape
    lon = float(icams.attrs['X_FIRST']) + float(icams.attrs['X_STEP'])*np.arange(sizex) - 0.5*float(icams.attrs['X_STEP'])
    lat = float(icams.attrs['Y_FIRST']) + float(icams.attrs['Y_STEP'])*np.arange(sizey) - 0.5*float(icams.attrs['Y_STEP'])
    icams = xr.DataArray(icams.tot_delay.values.reshape(sizey,sizex), coords=[lat, lon], dims=["lat", "lon"])
    ugeox = xr.DataArray(ugeo[0].values, coords=[ugeo.y.values, ugeo.x.values], dims=["lat", "lon"])
    icams = icams.interp_like(ugeox, method='linear')
    icams = icams.where(~np.isnan(ugeox))
    icams = icams.where(ugeox != 0) # just in case..
    icams = icams/ugeox  # converted to sltd [m]
    icams.values = mm2rad_s1(icams.values*1000)
    if outif:
        import lics_unwrap as lu
        lu.export_xr2tif(icams, outif, dogdal = False, set_to_pixel_registration = True)
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


# copied from lics_unwrap on 8th Oct 2025
def export_xr2tif(xrda, tif, lonlat=True, debug=True, dogdal=True, refto=None, set_to_pixel_registration=False):
    """Exports xarray dataarray to a geotiff

     Args:
        xrda (xarray.Dataarray): dataarray to export
        tif (string): path to output tif file
        lonlat (boolean): are the dimensions named as lon, lat?
        debug (boolean): just load it as float32
        dogdal (boolean): after exporting, perform gdalwarp (fix for potential issues in output geotiff) and gdal_translate to better compress
        refto (str): path to the (usually hgt) file to apply gdalwarp2match.py to. If None, it will not apply
        set_to_pixel_registration (boolean):  will rewrite header to assume Pixel Registration - that's by default in LiCSAR data (but not in rasterio...)
    """
    # coordsys = xrda.crs.split('=')[1]
    coordsys = "epsg:4326"
    if set_to_pixel_registration:
        dogdal = True  # pixel reg works correctly only through GDAL!
    if debug:
        xrda = xrda.astype(np.float32)
        # reset original spatial_ref
        if 'spatial_ref' in xrda:
            xrda = xrda.drop('spatial_ref')
        # remove attributes
        xrda.attrs = {}
    if lonlat:
        xrda = xrda.rename({'lon': 'x', 'lat': 'y'})
        # xrda = xrda.transpose('y', 'x')
    if xrda.y[1] > xrda.y[0]:
        xrda = xrda.sortby('y', ascending=False)
    xrda = xrda.rio.set_spatial_dims(x_dim="x", y_dim="y", inplace=True)
    #
    # Get pixel size
    dx = float(xrda['x'][1] - xrda['x'][0])
    dy = float(xrda['y'][0] - xrda['y'][1])
    # Get origin (top-left corner)
    x0 = xrda['x'][0].item() - dx / 2
    y0 = xrda['y'][0].item() - dy / 2
    # Build affine transform
    transform = Affine(dx, 0.0, x0, 0.0, -dy, y0)  # Note: dy is negative for north-up
    # Apply transform
    xrda.rio.write_transform(transform, inplace=True)
    #
    xrda = xrda.rio.write_crs(coordsys, inplace=True)
    if 'grid_mapping' in xrda.attrs:
        xrda.attrs.pop('grid_mapping', None)
    if dogdal:
        xrda.rio.to_raster(tif + 'tmp.tif')
        if refto:
            cmd = 'gdalwarp2match.py {0} {1} {2}; mv {2} {0}'.format(tif + 'tmp.tif', refto, tif)
            runcmd(cmd, printcmd=False)
        else:
            cmd = 'gdalwarp -t_srs EPSG:4326 {0} {1}'.format(tif + 'tmp.tif', tif)  # will fix some issues
            runcmd(cmd, printcmd=False)
        if set_to_pixel_registration:
            cmd = 'gdal_edit.py -mo AREA_OR_POINT=Point ' + tif
            runcmd(cmd, printcmd=False)
        cmd = 'mv {0} {1} 2>/dev/null; gdal_translate -of GTiff -co COMPRESS=DEFLATE -co PREDICTOR=3 {1} {0}'.format(
            tif, tif + 'tmp.tif')  # will compress better
        runcmd(cmd, printcmd=False)
        if os.path.exists(tif + 'tmp.tif'):
            os.remove(tif + 'tmp.tif')
        if os.path.exists(tif + 'tmp2.tif'):
            os.remove(tif + 'tmp2.tif')
    else:
        xrda.rio.to_raster(tif, compress='deflate')
        if set_to_pixel_registration:
            cmd = 'gdal_edit.py -mo AREA_OR_POINT=Point ' + tif
            runcmd(cmd, printcmd=False)


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




def fit_quadratic_xarray(da, degree = 'quadratic'):
    """ ready for either quadratic or cubic fitting"""
    # Get coordinates and flatten
    if 'lon' in da.coords:
        coordx='lon'
        coordy = 'lat'
    else:
        coordx = 'x'
        coordy = 'y'
    x, y = np.meshgrid(da[coordx], da[coordy])
    X = x.ravel()
    Y = y.ravel()
    Z = da.values.ravel()
    #
    # Mask NaNs
    mask = ~np.isnan(Z)
    X_valid, Y_valid, Z_valid = X[mask], Y[mask], Z[mask]
    #
    if degree == 'quadratic':
        # Design matrix for quadratic surface
        A = np.column_stack((X_valid**2, Y_valid**2, X_valid*Y_valid, X_valid, Y_valid, np.ones_like(X_valid)))
        model = sm.OLS(Z_valid, A).fit()
        #
        # Predict full surface
        A_full = np.column_stack((X**2, Y**2, X*Y, X, Y, np.ones_like(X)))
    elif degree == 'cubic':
        # Design matrix for full 2D cubic:
        # ax³ + by³ + cxy² + dx²y + ex² + fy² + gxy + hx + iy + j
        A = np.column_stack([
            X_valid ** 3, Y_valid ** 3,
            X_valid * Y_valid ** 2, X_valid ** 2 * Y_valid,
            X_valid ** 2, Y_valid ** 2,
            X_valid * Y_valid, X_valid, Y_valid,
            np.ones_like(X_valid)
        ])
        model = sm.OLS(Z_valid, A).fit()
        # Predict full surface
        A_full = np.column_stack([
            X ** 3, Y ** 3,
            X * Y ** 2, X ** 2 * Y,
            X ** 2, Y ** 2,
            X * Y, X, Y,
            np.ones_like(X)
        ])
    else:
        print('please use either quadratic or cubic degree')
        return
    # Predict full surface
    Z_fit = model.predict(A_full).reshape(da.shape)
    #
    return xr.DataArray(Z_fit, coords=da.coords, dims=da.dims)

import xarray as xr
import numpy as np
import statsmodels.api as sm

def fit_cubic_surface(da):
    # Flatten coordinates and values
    x, y = np.meshgrid(da.x, da.y)
    X = x.ravel()
    Y = y.ravel()
    Z = da.values.ravel()

    # Mask NaNs
    mask = ~np.isnan(Z)
    X_valid, Y_valid, Z_valid = X[mask], Y[mask], Z[mask]

    # Design matrix for full 2D cubic:
    # ax³ + by³ + cxy² + dx²y + ex² + fy² + gxy + hx + iy + j
    A = np.column_stack([
        X_valid**3, Y_valid**3,
        X_valid * Y_valid**2, X_valid**2 * Y_valid,
        X_valid**2, Y_valid**2,
        X_valid * Y_valid, X_valid, Y_valid,
        np.ones_like(X_valid)
    ])
    model = sm.OLS(Z_valid, A).fit()

    # Predict full surface
    A_full = np.column_stack([
        X**3, Y**3,
        X * Y**2, X**2 * Y,
        X**2, Y**2,
        X * Y, X, Y,
        np.ones_like(X)
    ])
    Z_fit = model.predict(A_full).reshape(da.shape)

    return xr.DataArray(Z_fit, coords=da.coords, dims=da.dims)


# based on ChatGPT5 - looks nice but... full of errors!
# from scipy.signal import convolve2d  # For the convolution

def calculate_mode(arr):
    """Calculates the mode of an array."""
    unique_values, counts = np.unique(arr, return_counts=True)
    return unique_values[np.argmax(counts)]


def mode_filter(data, window_size=(3, 3)):
    """Applies a mode filter using convolution to a 2D array.

    Args:
        data (np.ndarray): The 2D array to filter.
        window_size (tuple): The size of the mode filter window (rows, cols). MUST BE odd numbers, such as 3x3, 5x5,..

    Returns:
        np.ndarray: The filtered 2D array.
    """
    # Create a kernel (all ones) for the convolution
    # kernel = np.ones(window_size, dtype=np.int32)
    #
    # Pad the input array to handle edges.  This is VERY important to ensure
    # the output has the same dimensions as the input.  'same' padding adds
    # enough padding so the output is the same size as the input.  'reflect'
    # is usually a good choice for image/raster data.
    padded_data = np.pad(data, pad_width=((window_size[0] // 2, window_size[0] // 2),
                                          (window_size[1] // 2, window_size[1] // 2)),
                        mode='reflect')
    #
    # Perform the convolution.  Note: Using convolve2d directly can be slow for
    # very large arrays. For optimized performance, consider using libraries like
    # Dask or libraries specifically designed for image processing.
    # convolved = convolve2d(padded_data, kernel, mode='valid')  # 'valid' mode gives correct size after padding
    #
    # Calculate the mode for each window.
    mode_data = np.zeros_like(data)
    for i in range(data.shape[0]-int((window_size[0]-1)/2)):
        for j in range(data.shape[1] - int((window_size[0]-1)/2)):
            # Extract the window of data
            window = padded_data[i:i + window_size[0], j:j + window_size[1]]
            # Calculate and assign the mode
            mode_data[i+int((window_size[0]-1)/2), j+int((window_size[0]-1)/2)] = calculate_mode(window.flatten())
    #
    return mode_data



def apply_mode_filter_xr(raster, window_size=(3, 3)):
    """Applies a mode filter to a rioxarray DataArray.

    Args:
        raster (xr.DataArray): The rioxarray DataArray.
        window_size (tuple): The size of the mode filter window.

    Returns:
        xr.DataArray: The filtered rioxarray DataArray.
    """
    # Extract the data as a NumPy array
    data = raster.data[0]
    #
    # Apply the mode filter
    filtered_data = mode_filter(data, window_size)
    raster.values[0] = filtered_data
    return raster


def modefilter_mask(input_tiff, output_tiff):
    """Main function to read, filter, and write the GeoTIFF."""
    try:
        # Read the GeoTIFF using rioxarray
        raster = rioxarray.open_rasterio(input_tiff, masked=False)  # masked=False to handle 0/1 data

        # Apply the mode filter
        raster = apply_mode_filter_xr(raster, window_size=(3, 3))

        # Define compression options.  DEFLATE is a good general-purpose choice.
        compress_opts = {"compress": "DEFLATE", "zlevel": 5}  # zlevel controls compression level (1-9)

        # Write the filtered data to a new GeoTIFF with compression
        raster.rio.to_raster(output_tiff, **compress_opts)
        print(f"Successfully processed and saved to {output_tiff}")

    except Exception as e:
        print(f"An error occurred: {e}")