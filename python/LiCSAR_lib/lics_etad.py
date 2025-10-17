
# WORKING VERSION - means, it works.. but needs lot of improvements (that's first step only)
# e.g. interpolation can/must be improved

import s1etad
from s1etad import Sentinel1Etad, ECorrectionType

import os

import xarray as xr
import rioxarray
import geopandas as gpd
import numpy as np

from lics_unwrap import * # load_tif2xr
from shapely.ops import unary_union
from geocube.api.core import make_geocube
from scipy.interpolate import griddata

import matplotlib.pyplot as plt

def etacube_make_mask(eta, like):
    '''To create mask based on burst footprints in given grid.
    Args:
        eta (s1etad.Sentinel1Etad): to get the burst footprints from
        like (xr.Data*): datacube expected to have lon, lat coordinates
    '''
    bursts = []
    for aaa in eta.get_footprint():
        bursts.append(aaa)
    # merge burst polygons
    boundary = gpd.GeoSeries(unary_union(bursts))
    # and burn them as mask
    g = gpd.GeoDataFrame(
        {'mask' : [1]},
        index = [0],
        geometry=boundary,
        crs={"init": "epsg:4326"}
    )
    maskxr = make_geocube(vector_data=g, like=like.rio.set_spatial_dims(x_dim='lon',y_dim='lat'))
    maskxr=maskxr.rename({'x':'lon', 'y':'lat'})
    return maskxr


#from scipy.interpolate import LinearNDInterpolator
def interpolate2cube(merged_correction, like, method='linear'):
    '''Interpolation of the correction towards given data cube
    Args:
        merged_correction (dict): result of merged_correction function of s1etad. using only range (x)
        like (xr.Data*): datacube to read target dimensions from
        method (str): interpolation method ('nearest', 'linear', 'cubic')
    Note: using griddata triangulation - might be very slow indeed. but can be ok per burst?
    '''
    data = merged_correction['x']
    lat = merged_correction['lats']
    lon = merged_correction['lons']
    mask = data.mask
    values = data.data[~mask].ravel()
    lons = lon.data[~mask]
    lats = lat.data[~mask]
    #points=list(zip(lons, lats))
    points=list(zip(lons.ravel(), lats.ravel()))
    # initialize the linear interpolator
    #interp = LinearNDInterpolator(list(zip(lons, lats)), values)
    X, Y = np.meshgrid(like.lon.values, like.lat.values)
    #vals = interp(X, Y)
    vals = griddata(points, values, (X, Y), method=method)
    cube = like.copy()
    cube.values = vals
    return cube


def generate_correction_cube(eta, like, corrtype = 'ionospheric'):
    '''Function to generate datacube from the S1ETAD data
    Args:
        eta (s1etad.Sentinel1Etad): s1etad product - will use all bursts from here! Needs optimizations!
        like (xr.Data*): datacube to refer to (e.g. loaded geotiff to xr, to get dimensions etc.)
        corrtype (str): check help(ECorrectionType)
    '''
    print('extracting the correction from S1ETAD data')
    merged_correction = eta.merge_correction(corrtype, meter=True)
    print('interpolating '+corrtype+' correction')
    corrxr = interpolate2cube(merged_correction, like=like, method='linear')
    print('masking to burst extents')
    bmask = etacube_make_mask(eta, like=corrxr)
    corrxr = corrxr*bmask.mask
    return corrxr



def mm2rad_s1(inmm):
    speed_of_light = 299792458 #m/s
    radar_freq = 5.405e9  #for S1
    wavelength = speed_of_light/radar_freq #meter
    coef_r2m = -wavelength/4/np.pi*1000 #rad -> mm, positive is -LOS
    outrad = inmm/coef_r2m
    return outrad


def rad2mm_s1(inrad):
    speed_of_light = 299792458 #m/s
    radar_freq = 5.405e9  #for S1
    wavelength = speed_of_light/radar_freq #meter
    coef_r2m = -wavelength/4/np.pi*1000 #rad -> mm, positive is -LOS
    outmm = inrad*coef_r2m
    return outmm

etaf1 = '/gws/nopw/j04/nceo_geohazards_vol2/LiCS/temp/insar_proc/earmla/ETAD/etad/S1B_IW_ETA__AXDV_20180930T230604_20180930T230632_012950_017EB9_D6C7.SAFE'
etaf2 = '/gws/nopw/j04/nceo_geohazards_vol2/LiCS/temp/insar_proc/earmla/ETAD/etad/S1B_IW_ETA__AXDV_20181012T230605_20181012T230633_013125_018409_D9EC.SAFE'
eta1 = Sentinel1Etad(etaf1)
eta2 = Sentinel1Etad(etaf2)



bmask = etacube_make_mask(eta1, like=hgt)




#gacxr = (rad2mm_s1(gacoscorr)*bmask.mask*ifgxr.mask_extent)
#(-1*gacxr.where(gacxr != 0)).drop('spatial_ref').rename('mm').plot(vmin=-0, vmax=80)
#plt.title('GACOS correction: 20180930_20181012')

# full flow:
corrtype = 'ionospheric'
#corrtype = 'tropospheric'
#corrtype = 'sum' # not good!
#corrtype = 'geodetic' # perhaps this is SET?
#corrtype = 'doppler'

corrxr1 = generate_correction_cube(eta = eta1, like = hgt, corrtype = corrtype)
corrxr2 = generate_correction_cube(eta = eta2, like = hgt, corrtype = corrtype)
ifgc_m = corrxr2-corrxr1

#geodeticcorr=ifgc_m.copy()
#(1000*ifgc_m).drop('spatial_ref').rename('mm').plot()
#plt.title(corrtype+' correction: 20180930_20181012')
ionocorr=(ifgc_m*1000).copy()
#geodeticcorr=ifgc_m.copy()


# full flow:
corrtype = 'ionospheric'
corrtype = 'tropospheric'
#corrtype = 'sum' # not good!
#corrtype = 'geodetic' # perhaps this is SET?
#corrtype = 'doppler'

corrxr1 = generate_correction_cube(eta = eta1, like = hgt, corrtype = corrtype)
corrxr2 = generate_correction_cube(eta = eta2, like = hgt, corrtype = corrtype)
ifgc_m = corrxr2-corrxr1

#geodeticcorr=ifgc_m.copy()
#(1000*ifgc_m).drop('spatial_ref').rename('mm').plot()
#plt.title(corrtype+' correction: 20180930_20181012')
tropocorr=(ifgc_m*1000).copy()
#geodeticcorr=ifgc_m.copy()



# full flow:
corrtype = 'ionospheric'
corrtype = 'tropospheric'
corrtype = 'sum' # not good!
corrtype = 'geodetic' # perhaps this is SET?
#corrtype = 'doppler'
#hgt = ifgxr.hgt
corrxr1 = generate_correction_cube(eta = eta1, like = hgt, corrtype = corrtype)
corrxr2 = generate_correction_cube(eta = eta2, like = hgt, corrtype = corrtype)
ifgc_m = corrxr2-corrxr1

geodeticcorr=(ifgc_m*1000).copy()
(1000*ifgc_m).drop('spatial_ref').rename('mm').plot()
plt.title(corrtype+' correction: 20180930_20181012')

ifgxr=load_ifg(frame,pair,unw=False)
ifgxr

ifgxr['etad_solid']=geodeticcorr
ifgxr
...