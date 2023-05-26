import xarray as xr
from daz_iono import *
from lics_unwrap import *
from scipy.constants import speed_of_light
import numpy as np
from scipy.interpolate import griddata



def get_tecs_func(lat = 15.1, lon = 30.3, acq_times = [pd.Timestamp('2014-11-05 11:26:38'), pd.Timestamp('2014-11-29 11:26:38')]):
    return get_tecs(lat, lon, 800, acq_times[0:2], False)

def tecs_2d(xdf):
    func = lambda xdf: get_tecs(xdf.lat, xdf.lon, 800, acq_times[0:2], False)
    return xr.apply_ufunc(func, xdf)

def get_diff_tecs(lat = 15.1, lon = 30.3, acq_times = [pd.Timestamp('2014-11-05 11:26:38'), pd.Timestamp('2014-11-29 11:26:38')]):
    A,B = get_tecs(lat, lon, 800, acq_times[0:2], False)
    return B-A


def get_inc_frame(frame, heading=False):
    '''will get the incidence angle 2d xr.datarray
    if heading==True: return also heading raster
    '''
    metadir = os.path.join(os.environ['LiCSAR_public'], str(int(frame[:3])), frame, 'metadata')
    Ufile = os.path.join(metadir, frame + '.geo.U.tif')
    U = load_tif2xr(Ufile)
    U = U.where(U != 0)
    inc = U.copy()
    inc.values = np.degrees(np.arccos(U.values))
    if heading:
        Efile = os.path.join(metadir, frame + '.geo.E.tif')
        Nfile = os.path.join(metadir, frame + '.geo.N.tif')
        E = load_tif2xr(Efile)
        N = load_tif2xr(Nfile)
        E = E.where(E != 0)
        N = N.where(N != 0)
        heading = N.copy(deep=True)
        heading.values = np.degrees(np.arctan2(E, N)) + 90
        return inc, heading
    else:
        return inc


def get_resolution(ifg, in_m=True):
        """Gets resolution of the xr.dataset (or dataarray), either in metres or degrees
        """
        resdeg = (np.abs(hgt.lat[1]-hgt.lat[0])+np.abs(hgt.lon[1]-hgt.lon[0]))/2
        if in_m:
            latres = 111.32 * np.cos(np.radians(hgt.lat.mean())) * 1000 # in m
            return float(latres * resdeg)
        else:
            return float(resdeg)


# now for the test in licsar_disk/iono:
#ifgg=xr.open_dataset('ionotest.nc')
import matplotlib.pyplot as plt
#hgt = ifgg.hgt.where(ifgg.hgt != 0)




frame='149A_11032_131313'
frame='149A_11107_091009'
pair='20180930_20181012'


# load metadata, i.e.:

#master=20190504

#heading=-13.775063
#avg_incidence_angle=39.1918
#azimuth_resolution=14.068910
#range_resolution=2.329562
#avg_height=3165.176

epochs = pair.split('_')

# start using one epoch only
#acq = epochs[0]

# 1. get middle point - just super approx. for now
#hgt = ifg.hgt.where(ifg.hgt != 0)

#centre_range_m=880080.5691
#heading=-13.775063
#avg_incidence_angle=39.1918

inc=get_inc_frame(frame)
avg_incidence_angle = float(inc.mean())
# get hgt
metadir = os.path.join(os.environ['LiCSAR_public'],str(int(frame[:3])),frame,'metadata')
metafile = os.path.join(os.environ['LiCSAR_public'],str(int(frame[:3])),frame,'metadata','metadata.txt')
hgtfile=os.path.join(metadir, frame+'.geo.hgt.tif')
hgt = load_tif2xr(hgtfile)
hgt = hgt.where(hgt != 0)

scene_alt = float(hgt.median())
scene_center_lon = float(hgt.lon.mean())
scene_center_lat = float(hgt.lat.mean())

from LiCSAR_misc import *
center_time=grep1line('center_time',metafile).split('=')[1]
heading=float(grep1line('heading',metafile).split('=')[1])
centre_range_m=float(grep1line('centre_range_m',metafile).split('=')[1])
#center_time='23:06:29.844585'
#sat_alt_km = 800
#acqtime = pd.to_datetime(str(acq)+'T'+center_time)


def get_tecphase(epoch):
    acqtime = pd.to_datetime(str(epoch)+'T'+center_time)
    # this is to get point between sat and scene centre
    theta = np.radians(avg_incidence_angle)
    wgs84 = nv.FrameE(name='WGS84')
    Pscene_center = wgs84.GeoPoint(latitude=scene_center_lat, longitude=scene_center_lon, degrees=True)
    #    burst_len = 7100*2.758277 #approx. satellite velocity on the ground 7100 [m/s] * burst_interval [s]
        ###### do the satg_lat, lon
    azimuthDeg = heading-90 #yes, azimuth is w.r.t. N (positive to E)
    elevationDeg = 90-avg_incidence_angle # this is to get the avg sat altitude/range
    slantRange = centre_range_m
    # from daz_iono:
    x, y, z = aer2ecef(azimuthDeg, elevationDeg, slantRange, scene_center_lat, scene_center_lon, scene_alt)
    satg_lat, satg_lon, sat_alt = ecef2latlonhei(x, y, z)
    sat_alt_km = round(sat_alt/1000)
    Psatg = wgs84.GeoPoint(latitude=satg_lat, longitude=satg_lon, degrees=True)
    # get middle point between scene and sat - and get F2 height for it
    path = nv.GeoPath(Pscene_center.to_nvector(), Psatg.to_nvector())
    # get point in the middle
    Pmid_scene_sat = path.interpolate(0.5).to_geo_point()
    # get hionos in that middle point:
    tecs, hionos = get_tecs(Pmid_scene_sat.latitude_deg, Pmid_scene_sat.longitude_deg, sat_alt_km, [acqtime], returnhei = True)
    hiono = hionos[0]*1000 # m
    # first, get IPP - ionosphere pierce point
    # range to IPP can be calculated using:
    range_IPP = slantRange * hiono / sat_alt
    #
    # so now let's get the IPP coordinates, using the range to IPP --- BUT, first we need to update the elevationDeg, as the 
    # ionospheric plasma would have similar effect to the projected scene as your leg projected inside water w.r.t. outside (a 'cut' appears, i.e. change in look angle)
    # get inc angle at IPP - see iono. single layer model function
    #earth_radius = 6378160 # m
    #sin_thetaiono = earth_radius/(earth_radius+hiono) * np.sin(theta)
    x, y, z = aer2ecef(azimuthDeg, elevationDeg, range_IPP, scene_center_lat, scene_center_lon, scene_alt)
    ippg_lat, ippg_lon, ipp_alt = ecef2latlonhei(x, y, z)
    #
    dlat = ippg_lat-scene_center_lat
    dlon = ippg_lon-scene_center_lon
    #
    # now i need to shift all the points towards the satellite, by the path_scenecenter_to_IPP distance (direction)
    #
    # this is to grid to less points:
    ionosampling=10000 # m  --- by default, 10 km sampling should be ok?
    resolution = get_resolution(hgt, in_m=True)  # just mean avg in both lon, lat should be ok
    # how large area is covered
    lonextent = len(hgt.lon)*resolution
    # so what is the multilook factor?
    mlfactorlon = round(len(hgt.lon)/(lonextent/ionosampling))
    latextent = len(hgt.lat)*resolution
    mlfactorlat = round(len(hgt.lat)/(latextent/ionosampling))
    hgtml = hgt.coarsen({'lat': mlfactorlat, 'lon': mlfactorlon}, boundary='trim').mean()
    incml = inc.coarsen({'lat': mlfactorlat, 'lon': mlfactorlon}, boundary='trim').mean()
    # get range towards iono single-layer in the path to the satellite
    range2iono = (hiono - hgtml) / np.cos(np.radians(incml))
    earth_radius = 6378160  # m
    ionoxr = incml.copy(deep=True)
    print('getting TEC values sampled by {} km. in latitude:'.format(str(round(ionosampling / 1000))))
    for i in range(len(range2iono.lat.values)):
        print(str(i) + '/' + str(len(range2iono.lat.values)))
        for j in range(len(range2iono.lon.values)):
            if ~np.isnan(incml.values[i, j]):
                #theta = float(np.radians(incml.values[i, j]))
                eledeg = float(90 - incml.values[i, j])
                ilat_ground, ilon_ground = range2iono.lat.values[i], range2iono.lon.values[j]
                x, y, z = aer2ecef(azimuthDeg, eledeg, range2iono.values[i, j], ilat_ground, ilon_ground, float(hgtml.values[i,j]))
                ilat, ilon, ialt = ecef2latlonhei(x, y, z)
                theta = float(np.radians(incml.values[i, j]))
                sin_thetaiono = earth_radius / (earth_radius + hiono) * np.sin(theta)
                ionoxr.values[i, j] = get_tecs(ilat, ilon, sat_alt_km, [acqtime], False)[0] / np.sqrt(1 - sin_thetaiono ** 2) # maybe no need for the last term??
    '''
    # so now use the shifted coordinates of decimated hgt to find tec:
    #fortec_lon = incml.lon.values + dlon
    #fortec_lat = incml.lat.values + dlat
    # actually, we can do it (much) better than that! see above... test
    #
    #
    # do it the good old way:
    # also correct for geometric squinting through ionosphere, at the Hiono:
    earth_radius = 6378160  # m
    ionoxr = incml.copy(deep=True)
    print('getting TEC values sampled by {} km. in latitude:'.format(str(round(ionosampling/1000))))
    for i in range(len(fortec_lat)):
        print(str(i)+'/'+str(len(fortec_lat)))
        for j in range(len(fortec_lon)):
            ilat, ilon = fortec_lat[i], fortec_lon[j]
            theta = float(np.radians(incml.values[i,j]))
            # now for the iono-squint fix:
            sin_thetaiono = earth_radius / (earth_radius + hiono) * np.sin(theta)
            ionoxr.values[i,j] = get_tecs(ilat, ilon, sat_alt_km, [acqtime], False)[0]/np.sqrt(1-sin_thetaiono**2)
    #
    '''
    # ok, now correct for geometric squinting through ionosphere, at the Hiono:
    #sin_thetaiono = earth_radius/(earth_radius+hiono) * np.sin(theta)
    #ionoxr = ionoxr/np.sqrt(1-sin_thetaiono**2)
    #def mm2rad_s1(inmm):
    #speed_of_light = 299792458 #m/s
    #radar_freq = 5.405e9  #for S1
    #wavelength = speed_of_light/radar_freq #meter
    #coef_r2m = -wavelength/4/np.pi*1000 #rad -> mm, positive is -LOS
    #outrad = inmm/coef_r2m
    #return outrad
    # now, convert TEC values into 'phase' - simplified here (?)
    f0 = 5.4050005e9
    #inc = avg_incidence_angle  # e.g. 39.1918 ... oh but... it actually should be the iono-squint-corrected angle. ignoring now
    ionoxr = -4*np.pi*40.308193/speed_of_light/f0*ionoxr/np.cos(np.radians(incml))
    #ionoxr = -2 * np.pi * 40.308193 / speed_of_light / f0 * ionoxr / np.cos(np.radians(incml))
    # ionodelay in seconds: dT=2*40.308193/speed_of_light/f0^2 * ionoxr/np.cos(np.radians(incml))
    # so the diff phase would be: pha[rad] = -2 pi * dT * f0;
    # while distance would be: d[m] = pha * -speed_of_light/f0 /(4 pi) = 0.5 dT * speed_of_light       ; v=s/t -> d=c*dT <- in both ways
    # BUT! i might have wrong IPP angle values - just because i don't scale lons, only shift them!
    # now the ionoxr contains phase in radians
    return ionoxr

# get tec phase for both epochs:
tecphase1=get_tecphase(epochs[0])
tecphase2=get_tecphase(epochs[1])
# and their difference
tecdiff = tecphase2-tecphase1

# so the final step is to interpolate the result - use bilinear int. - and get it back to hgt (or ifg, perform correction)
tecout=tecdiff.interp(lat=inc.lat.values, lon=inc.lon.values, method="linear",kwargs={"fill_value": "extrapolate"})
tecout.to_netcdf('tecout.nc')

# 2022-10-17:
# something is missing - maybe the effect of range offsets due to ionosphere, mapped by range px offset tracking during coreg towards M?
master=str(grep1line('master',metafile).split('=')[1])
tecphaseM=get_tecphase(m)
tecdiff2 = tecphaseM-tecphase1
tecout2=tecdiff2.interp(lat=inc.lat.values, lon=inc.lon.values, method="linear",kwargs={"fill_value": "extrapolate"})
tecdiff3 = tecphaseM-tecphase2
tecout3=tecdiff3.interp(lat=inc.lat.values, lon=inc.lon.values, method="linear",kwargs={"fill_value": "extrapolate"})


#frame='069A_06694_131402'
#pair='20170304_20170316'

ifg=load_ifg(frame,pair)
bagr=ifg.pha.copy()
bagr.values = wrap2phase(ifg.pha.values - tecout.values - tecout2.values - tecout3.values)
bagr.plot()

# done!!!

'''

def interpolate_tecdiff2ifg(tecdiff, ifgxr): 
    '''ifgxr must be just xr.dataarray!
    '''
    teclats = tecdiff.lat.values
    teclons = tecdiff.lon.values
    values = tecdiff.values.ravel()
    points=list(zip(teclons, teclats))
    X, Y = np.meshgrid(ifgxr.lon.values, ifgxr.lat.values)
    vals = griddata(points, values, (X, Y), method='linear')
    cube = ifgxr.copy()
    cube.values = vals
    return cube




dsi = ds.interp(lat=new_lat, lon=new_lon)


ifgg['ionocorr']=interpolate_tecdiff2ifg(tecdiff, ifgxr)





data = merged_correction['x']
    lat = merged_correction['lats']
    lon = merged_correction['lons']
    mask = data.mask
    values = data.data[~mask].ravel()
    lons = lon.data[~mask]
    lats = lat.data[~mask]
    points=list(zip(teclons, teclats))

    # initialize the linear interpolator
    #interp = LinearNDInterpolator(list(zip(lons, lats)), values)
    X, Y = np.meshgrid(like.lon.values, like.lat.values)
    #vals = interp(X, Y)
    vals = griddata(points, values, (X, Y), method=method)
    cube = like.copy()
    cube.values = vals
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
    
    
    

a['tecphase'] = -4*np.pi*40.308193/speed_of_light/f0*a['tecdiffs']/np.cos(np.radians(inc))
a['tecphase'].plot()





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
'''