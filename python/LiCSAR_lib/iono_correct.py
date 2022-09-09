import xarray as xr
from daz_iono import *


def get_tecs_func(lat = 15.1, lon = 30.3, acq_times = [pd.Timestamp('2014-11-05 11:26:38'), pd.Timestamp('2014-11-29 11:26:38')]):
    return get_tecs(lat, lon, 800, acq_times[0:2], False)

def tecs_2d(xdf):
    func = lambda xdf: get_tecs(xdf.lat, xdf.lon, 800, acq_times[0:2], False)
    return xr.apply_ufunc(func, xdf)

def get_diff_tecs(lat = 15.1, lon = 30.3, acq_times = [pd.Timestamp('2014-11-05 11:26:38'), pd.Timestamp('2014-11-29 11:26:38')]):
    A,B = get_tecs(lat, lon, 800, acq_times[0:2], False)
    return B-A



# load metadata, i.e.:

master=20190504

heading=-13.775063
avg_incidence_angle=39.1918
azimuth_resolution=14.068910
range_resolution=2.329562

avg_height=3165.176


frame='149A_11032_131313'
pair='20180930_20181012'

epochs = pair.split('_')
acq = epochs[0]

# 1. get middle point - just super approx. for now
hgt = ifg.hgt.where(ifg.hgt != 0)
scene_alt = float(hgt.median())
scene_center_lon = float(hgt.lon.mean())
scene_center_lat = float(hgt.lat.mean())
centre_range_m=880080.5691
heading=-13.775063
avg_incidence_angle=39.1918


sat_alt_km = 800
acqtime = pd.to_datetime(str(acq)+'T'+center_time)

# this is to get point between sat and scene centre
theta = np.radians(avg_incidence_angle)
wgs84 = nv.FrameE(name='WGS84')
Pscene_center = wgs84.GeoPoint(latitude=scene_center_lat, longitude=scene_center_lon, degrees=True)
#    burst_len = 7100*2.758277 #approx. satellite velocity on the ground 7100 [m/s] * burst_interval [s]
    ###### do the satg_lat, lon
azimuthDeg = heading-90 #yes, azimuth is w.r.t. N (positive to E)
elevationDeg = 90-avg_incidence_angle
slantRange = centre_range_m
# from daz_iono:
x, y, z = aer2ecef(azimuthDeg, elevationDeg, slantRange, scene_center_lat, scene_center_lon, scene_alt)
satg_lat, satg_lon, sat_alt = ecef2latlonhei(x, y, z)
Psatg = wgs84.GeoPoint(latitude=satg_lat, longitude=satg_lon, degrees=True)
# get middle point between scene and sat - and get F2 height for it
path = nv.GeoPath(Pscene_center.to_nvector(), Psatg.to_nvector())
# get point in the middle
Pmid_scene_sat = path.interpolate(0.5).to_geo_point()
# get hionos in that middle point:
tecs, hionos = get_tecs(Pmid_scene_sat.latitude_deg, Pmid_scene_sat.longitude_deg, sat_alt_km, [acqtime], returnhei = True)

########################## CHECK THIS BELOW:
hiono = hionos[0]*1000 # m
tec = tecs[0]

# first, get IPP - ionosphere pierce point
# range to IPP can be calculated using:
range_IPP = hiono/np.sin(theta)
x, y, z = aer2ecef(azimuthDeg, elevationDeg, range_IPP, scene_center_lat, scene_center_lon, scene_alt)
ippg_lat, ippg_lon, ipp_alt = ecef2latlonhei(x, y, z)
Pippg = wgs84.GeoPoint(latitude=ippg_lat, longitude=ippg_lon, degrees=True)
path_scenecenter_Pippg = nv.GeoPath(Pscene_center, Pippg)

# now i need to shift all the points towards the satellite, by the path_scenecenter_to_IPP distance (direction)
neco jako... displace..with geopath...
def displace_scene(scene, geopath):
    return newscene


ionoscene = hgt.copy()
ionoscene = displace_scene(ionoscene)






def get_midpoint_towards_satellite(lon, lat, gpoint):
    gpoint2 = wgs84.GeoPoint(latitude=lat, longitude=lon, degrees=True)
    # get middle point between scene and sat - and get F2 height for it
    path = nv.GeoPath(gpoint2.to_nvector(), gpoint.to_nvector())
    # get point in the middle - might be not that accurate, but perhaps ok for range direction...
    midpoint = path.interpolate(0.5).to_geo_point()
    return midpoint



TECV_B = get_tecs(PippB.latitude_deg, PippB.longitude_deg, round(sat_alt/1000), [epochdate], False)[0]
    # get inc angle at IPP - see iono. single layer model function
    earth_radius = 6378160 # m
    sin_thetaiono = earth_radius/(earth_radius+hiono) * np.sin(theta)
    TECS_A = TECV_A/np.sqrt(1-sin_thetaiono**2)
#








# this is to grid to less points:
ionosampling=5000 # m
from lics_unwrap import *
resolution = get_resolution(hgt, in_m=True)
# how large area is covered
lonextent = len(hgt.lon)*resolution
# so how many pixels do we want?
mlfactor = round(len(hgt.lon)/(lonextent/ionosampling))
hgtml = hgt.coarsen({'lat': mlfactor, 'lon': mlfactor}, boundary='trim').mean()



# this is to get TECs for the acquisition times:
sat_alt_km = 800
epochs = pair.split('_')
for acq in epochs:
    acqtime = pd.to_datetime(str(acq)+'T'+center_time)
    get_tecs(Pmid_scene_sat.latitude_deg, Pmid_scene_sat.longitude_deg, sat_alt_km, [acqtime], returnhei = False)
center_time='23:06:29.844585'
# get hionos in that middle point:
tecs, hionos = get_tecs(Pmid_scene_sat.latitude_deg, Pmid_scene_sat.longitude_deg, 800, acq_times, returnhei = True)




i=0
a['tecdiffs'] = a.z.copy()
for ilat in range(len(a.lat)):
    glat = a.lat.values[ilat]
    i=i+1
    print('done lat '+str(i)+' from '+str(len(a.lat.values)))
    for ilon in range(len(a.lon)):
        glon = a.lon.values[ilon]
        tdiff = tecs_diff(glat, glon)
        a['tecdiffs'].values[ilat,ilon] = tdiff



# last - get the phase correction!


from scipy.constants import speed_of_light
f0 = 5.4050005e9
inc = avg_incidence_angle  # e.g. 39.1918
a['tecphase'] = -4*np.pi*40.308193/speed_of_light/f0*a['tecdiffs']/np.cos(np.radians(inc))
a['tecphase'].plot()


# and finally interpolate it to the whole grid (use griddata - linear?)


a = hgt.sel(lat=slice(None,None,100),lon=slice(None,None,100))


b = xr.apply_ufunc(func, a.lat, a.lon)


ifg = xr.open_dataset('ionotest.nc')






i=0
a['tecdiffs'] = a.z.copy()
for ilat in range(len(a.lat)):
    glat = a.lat.values[ilat]
    i=i+1
    print('done lat '+str(i)+' from '+str(len(a.lat.values)))
    for ilon in range(len(a.lon)):
        glon = a.lon.values[ilon]
        tdiff = tecs_diff(glat, glon)
        a['tecdiffs'].values[ilat,ilon] = tdiff






center_time=23:06:29.844585
heading=-13.775063


hgt = hgt.rename({'x':'lon', 'y':'lat'})
a = hgt.sel(lat=slice(None,None,100),lon=slice(None,None,100))







frame = '099A_05417_131313'
hgtif = '/home/home02/earmla/licsar/eq/frames/099A_05417_131313/GEOC/geo/099A_05417_131313.geo.hgt.tif'

master='2016-09-07'
dtime='11:26:38'
heading=-10.140217
sathei = 728 #km


import rioxarray
import iri2016 # as iri
import numpy as np
import pandas as pd
import os
from scipy.constants import speed_of_light
from lics_unwrap import load_tif2xr


hgttif = '/gws/nopw/j04/nceo_geohazards_vol1/public/LiCSAR_products/149/149A_11032_131313/metadata/149A_11032_131313.geo.hgt.tif'
hgt = load_tif2xr(hgttif)




master=20190504
center_time=23:06:29.844585
heading=-13.775063
avg_incidence_angle=39.1918
azimuth_resolution=14.068910
range_resolution=2.329562
centre_range_m=880080.5691
avg_height=3165.176






outfile='test_gacos_ifg.tif'
#make_gacos_ifg(frame, pair, outfile)
gacoscorr = load_tif2xr(outfile)



eptxt = '/home/home02/earmla/licsar/eq/frames/099A_05416_131313/epochs.txt'
epochs = pd.read_csv(eptxt, header=None)
acq_times = []
for e in epochs[0]:
    acq_times.append(pd.to_datetime(str(e)+'T'+dtime))
#acq_times
