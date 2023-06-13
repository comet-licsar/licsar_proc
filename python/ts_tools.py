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


def stamps2cum(datesfile='date.txt', llfile='ps_ll.txt', uwphasefile='uw_phase', velsfile='vels.txt', wavelength=0.0555, resol=0.00075, outncfile = 'ers_cumm.nc'):
    """ Converts standard stamps outputs to the 'cum' table.
    """
    print('loading data')
    dates=pd.read_csv(datesfile, header=None)
    # convert dates to %y-%m-%d
    dates=dates.apply(lambda x: str(pd.Timestamp(str(int(x))).date()), axis=1)
    lonlat=pd.read_csv(llfile, header=None, delimiter='   ', engine='python')
    lonlat=lonlat.rename(columns={0: "LON", 1: "LAT"})
    uwpha=pd.read_csv(uwphasefile, header=None, delimiter='  ', engine='python')
    uwpha=pd.DataFrame(uwpha.values, columns=dates.values)
    # convert to mm, as expected in 'cum':
    coef_r2m = -wavelength/4/np.pi*1000 #rad -> mm
    cum = coef_r2m * uwpha
    uwpha = ''
    cum[['LON','LAT']]=lonlat
    if velsfile:
        vels=pd.read_csv(velsfile, header=None, delimiter='  ', engine='python')[2]
        print('not ready yet, i did it manually this way:')
        '''
        df=df.rename(columns={2:'vel'})
        df[['LON','LAT']]=lonlat
        ... stuff below... and then:
        anc['vel']=nc
        anc.rio.set_spatial_dims(x_dim="lon", y_dim="lat", inplace=True)
        anc.rio.set_spatial_dims(x_dim="lon", y_dim="lat", inplace=True)
        encode = {'vel': {'zlib': True, 'complevel': 9},
 'cum': {'zlib': True, 'complevel': 9},
           'time': {'dtype': 'i4'}
          }

        anc.to_netcdf(outncfile, encoding=encode)
        '''
    # from below:
    print('converting to grid')
    df=cum
    to_bin = lambda x: np.floor(x / resol) * resol
    df["lat"] = to_bin(df['LAT'])
    df["lon"] = to_bin(df['LON'])
    groups = df.groupby(["lat", "lon"])
    medgrid = groups.agg(np.nanmedian)
    lat = medgrid.index.get_level_values(level=0)
    lon = medgrid.index.get_level_values(level=1)
    dates=dates.to_list()
    dates.sort()
    
    # now convert from dates
    datum = dates[0]
    a=medgrid[datum].to_xarray().assign_coords(
        {'time':dt.datetime.strptime(datum, '%Y-%m-%d')}).expand_dims('time').rename('cum')
    for datum in dates[1:]:
        b = medgrid[datum].to_xarray().assign_coords(
            {'time':dt.datetime.strptime(datum, '%Y-%m-%d')}).expand_dims('time').rename('cum')
        a = xr.concat([a,b], dim='time')
    a = a - a[0] # make the first epoch zero
    
    print('storing')
    nc=xr.Dataset()
    nc['cum']=a
    #nc = nc.assign_coords({'time':a.time.values})
    # set ref. system
    nc=nc.rio.set_spatial_dims(x_dim="lon", y_dim="lat", inplace=True)
    nc=nc.rio.write_crs("EPSG:4326", inplace=True)
    encode = {
          'cum': {'zlib': True, 'complevel': 9},
          'time': {'dtype': 'i4'}
         }
    if outncfile:
        nc.to_netcdf(outncfile, encoding=encode)
    return nc


# resol_m is approx resol*111111 [m]
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
