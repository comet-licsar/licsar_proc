#!/usr/bin/env python3
import pandas as pd
import numpy as np
import xarray as xr
import rioxarray
import datetime as dt
import matplotlib.pyplot as plt

# 0.000274 for 30 m
def csv2nc(csv, outncfile = None, resol = 0.0025, extracols = []):
    """ Converts a csv file into netcdf.
    The CSV must have following columns with the header as the first line like as follows (order does not matter, but keep the upper case):
    LAT,LON,VEL,COHER,*dates*
    where *dates* are column names in format of yyyy-mm-dd.
    
    Usage example:
        nc = csv2nc('test.csv', 'test.nc')
        nc.vel.plot(); plt.show()
    
    Args:
        csv (string): path to the csv
        outncfile (string): path to the output netcdf file (or None to only return the object)
        resol (float): output resolution cell size in degrees (WGS-84)
        extracols (list): list of column names to be also exported (apart from VEL,COHER and dates)
    Returns:
        string: path to generated snaphu.conf
    """
    df = pd.read_csv(csv)
    print('converting to the netcdf file')
    nc = df2nc(df, outncfile = outncfile, resol = resol, extracols = extracols)
    return nc


def df2nc(df, outncfile = None, resol = 0.0025, extracols = []): #, cumthecum = False):
    """ Converts pandas dataframe (loaded csv file) to NetCDF.
    See help of csv2nc for proper formatting.
    Grid would aggregate values in each cell by their median.
    """
    to_bin = lambda x: np.floor(x / resol) * resol
    df["lat"] = to_bin(df['LAT'])
    df["lon"] = to_bin(df['LON'])
    groups = df.groupby(["lat", "lon"])
    medgrid = groups.agg(np.nanmedian)
    
    lat = medgrid.index.get_level_values(level=0)
    lon = medgrid.index.get_level_values(level=1)
    dates = df.columns[df.columns.str.match(r"\d{4}-\d{2}-\d{2}")].to_list()
    dates.sort()
    
    #cols = ['VEL_U', 'COHER'] #, 'SIGMA VEL_U']
    cols = ['VEL', 'COHER']
    cols = cols + extracols
    nc = medgrid[cols].to_xarray()
    #velcol = 'vel'
    nc = nc.rename({'VEL': 'vel','COHER': 'coh'})
    
    # now convert from dates
    datum = dates[0]
    a=medgrid[datum].to_xarray().assign_coords(
        {'time':dt.datetime.strptime(datum, '%Y-%m-%d')}).expand_dims('time').rename('cum')
    for datum in dates[1:]:
        b = medgrid[datum].to_xarray().assign_coords(
            {'time':dt.datetime.strptime(datum, '%Y-%m-%d')}).expand_dims('time').rename('cum')
        a = xr.concat([a,b], dim='time')
    a = a - a[0] # make the first epoch zero
    #if cumthecum:
    #    a = a.cumsum()
    nc = nc.assign_coords({'time':a.time.values})
    nc['cum'] = a
    # set ref. system
    nc.rio.set_spatial_dims(x_dim="lon", y_dim="lat", inplace=True)
    nc.rio.write_crs("EPSG:4326", inplace=True)
    #compress it and store as netcdf
    encode = {'vel': {'zlib': True, 'complevel': 9},
          'coh': {'zlib': True, 'complevel': 9},
          'cum': {'zlib': True, 'complevel': 9},
          'time': {'dtype': 'i4'}
         }
    if outncfile:
        nc.to_netcdf(outncfile, encoding=encode)
    return nc


'''
def invert_ts_custom(cube, timestamps):
    
'''
