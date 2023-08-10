import xarray as xr
from daz_iono import *
from lics_unwrap import *
from scipy.constants import speed_of_light
import numpy as np
from scipy.interpolate import griddata
from LiCSAR_misc import *


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


def get_resolution(hgt, in_m=True):
        """Gets resolution of the xr.dataset (or dataarray), either in metres or degrees
        """
        resdeg = (np.abs(hgt.lat[1]-hgt.lat[0])+np.abs(hgt.lon[1]-hgt.lon[0]))/2
        if in_m:
            latres = 111.32 * np.cos(np.radians(hgt.lat.mean())) * 1000 # in m
            return float(latres * resdeg)
        else:
            return float(resdeg)



def make_ionocorr_pair(frame, pair, source = 'code', outif=None):
    """ This will generate ionospheric correction for given frame-pair.
    It would optionally output the result to a geotiff.
    
    Args:
        frame (str):    frame ID
        pair (str):     pair (e.g. '20180930_20181012')
        source (str):   source model for TEC values. Either 'iri' or 'code'.
        outif (str):    if given, will export the iono phase screen to given geotiff
    Returns:
        xr.DataArray:   estimated ionospheric phase screen
    """
    ifg = load_ifg(frame, pair)
    epochs = pair.split('_')
    #
    tecphase1 = make_ionocorr_epoch(frame, epochs[0], source = source)
    tecphase2 = make_ionocorr_epoch(frame, epochs[1], source = source)
        # and their difference
    tecdiff = tecphase1 - tecphase2  # 07/2023: should it be this way/opposite???!!!! (i think so)
    #    # tecdiff = interpolate_nans_pyinterp(tecdiff)
    #tecdiff = interpolate_nans_bivariate(tecdiff)
    #tecdiff = tecdiff.interp_like(ifg, method='linear', kwargs={"bounds_error": False, "fill_value": None})
    #tecdiff = interpolate_nans_bivariate(tecdiff) # not needed?
    #    if np.max(np.isnan(tecdiff.values)):
    #        tecdiff = interpolate_nans_bivariate(tecdiff)
    tecdiff = tecdiff.where(ifg.mask_extent == 1)
    if outif:
        export_xr2tif(tecdiff,outif)
    return tecdiff



def make_ionocorr_epoch(frame, epoch, source = 'code', fixed_f2_height_km = None):
    """
    Args:
        ...
        fixed_f2_height_km: if None, it will use IRI2016 to estimate the height (in mid point between scene centre and satellite)
    """
    #if source == 'code':
    #    # this is to grid to less points:
    #    ionosampling=10000 # m 
    #else:
    ionosampling = 20000  # m  --- by default, 20 km sampling should be ok?
    #if source == 'code':
    #    ionosampling=20000 # m  --- by default, 20 km sampling should be ok?
    #else:
    #    ionosampling=40000
    # start using one epoch only
    #acq = epochs[0]
    #
    # 1. get middle point - just super approx. for now
    #hgt = ifg.hgt.where(ifg.hgt != 0)
    #
    #centre_range_m=880080.5691
    #heading=-13.775063
    #avg_incidence_angle=39.1918
    #
    inc=get_inc_frame(frame)
    avg_incidence_angle = float(inc.mean())
    # get hgt # no need anymore (2023/08)
    metadir = os.path.join(os.environ['LiCSAR_public'],str(int(frame[:3])),frame,'metadata')
    metafile = os.path.join(os.environ['LiCSAR_public'],str(int(frame[:3])),frame,'metadata','metadata.txt')
    #hgtfile=os.path.join(metadir, frame+'.geo.hgt.tif')
    #hgt = load_tif2xr(hgtfile)
    #hgt = hgt.where(hgt != 0)
    #
    #scene_alt = float(hgt.median())
    scene_center_lon = float(inc.lon.mean())
    scene_center_lat = float(inc.lat.mean())
    #
    center_time=grep1line('center_time',metafile).split('=')[1]
    heading=float(grep1line('heading',metafile).split('=')[1])
    centre_range_m=float(grep1line('centre_range_m',metafile).split('=')[1])
    #
    master=str(grep1line('master',metafile).split('=')[1])
    #
    acqtime = pd.to_datetime(str(epoch) + 'T' + center_time)
    # this is to get point between sat and scene centre
    theta = np.radians(avg_incidence_angle)
    wgs84 = nv.FrameE(name='WGS84')
    Pscene_center = wgs84.GeoPoint(latitude=scene_center_lat, longitude=scene_center_lon, degrees=True)
    #    burst_len = 7100*2.758277 #approx. satellite velocity on the ground 7100 [m/s] * burst_interval [s]
    ###### do the satg_lat, lon
    azimuthDeg = heading - 90  # yes, azimuth is w.r.t. N (positive to E)
    elevationDeg = 90 - avg_incidence_angle  # this is to get the avg sat altitude/range
    slantRange = centre_range_m
    # from daz_iono:
    x, y, z = aer2ecef(azimuthDeg, elevationDeg, slantRange, scene_center_lat, scene_center_lon, 0) #scene_alt) ### this is wrt ellipsoid!
    satg_lat, satg_lon, sat_alt = ecef2latlonhei(x, y, z)
    sat_alt_km = round(sat_alt / 1000)
    Psatg = wgs84.GeoPoint(latitude=satg_lat, longitude=satg_lon, degrees=True)
    # get middle point between scene and sat - and get F2 height for it
    path = nv.GeoPath(Pscene_center.to_nvector(), Psatg.to_nvector())
    # get point in the middle
    Pmid_scene_sat = path.interpolate(0.5).to_geo_point()
    if fixed_f2_height_km:
        #hiono = 450
        hiono = fixed_f2_height_km
    else:
        # get estimated hiono from IRI2016 in that middle point (F2 peak altitude):
        try:
            tecs, hionos = get_tecs(Pmid_scene_sat.latitude_deg, Pmid_scene_sat.longitude_deg, sat_alt_km, [acqtime],
                            returnhei=True)
            hiono = hionos[0]
        except:
            print('error in IRI2016, perhaps not installed? Setting standard F2 peak altitude')
            hiono = 450
    print('Getting IPP in the altitude of {} km'.format(str(int(hiono))))
    hiono = hiono * 1000  # m
    # first, get IPP - ionosphere pierce point
    # range to IPP can be calculated using:
    range_IPP = slantRange * hiono / sat_alt
    #
    # so now let's get the IPP coordinates, using the range to IPP --- BUT, first we need to update the elevationDeg, as the
    # ionospheric plasma would have similar effect to the projected scene as your leg projected inside water w.r.t. outside (a 'cut' appears, i.e. change in look angle)
    # get inc angle at IPP - see iono. single layer model function
    # earth_radius = 6378160 # m
    # sin_thetaiono = earth_radius/(earth_radius+hiono) * np.sin(theta)
    x, y, z = aer2ecef(azimuthDeg, elevationDeg, range_IPP, scene_center_lat, scene_center_lon, 0) #scene_alt) # this wrt ellipsoid
    ippg_lat, ippg_lon, ipp_alt = ecef2latlonhei(x, y, z)
    #
    dlat = ippg_lat - scene_center_lat
    dlon = ippg_lon - scene_center_lon
    #
    # now i need to shift all the points towards the satellite, by the path_scenecenter_to_IPP distance (direction)
    #
    resolution = get_resolution(inc, in_m=True)  # just mean avg in both lon, lat should be ok
    # how large area is covered
    lonextent = len(inc.lon) * resolution
    # so what is the multilook factor?
    mlfactorlon = round(len(inc.lon) / (lonextent / ionosampling))
    latextent = len(inc.lat) * resolution
    mlfactorlat = round(len(inc.lat) / (latextent / ionosampling))
    #hgtml = hgt.coarsen({'lat': mlfactorlat, 'lon': mlfactorlon}, boundary='trim').mean()
    incml = inc.coarsen({'lat': mlfactorlat, 'lon': mlfactorlon}, boundary='trim').mean()
    # get range towards iono single-layer in the path to the satellite, consider hgt
    # range2iono = (hiono - hgtml) / np.cos(np.radians(incml))
    # get range towards iono single-layer in the path to the satellite, do not consider hgt
    range2iono = hiono / np.cos(np.radians(incml))
    earth_radius = 6378160  # m
    ionoxr = incml.copy(deep=True)
    if source == 'code':
        tecxr=get_vtec_from_code(acqtime, lat=0, lon=0, return_fullxr = True)
    print('getting TEC values sampled by {} km.'.format(str(round(ionosampling / 1000))))
    for i in range(len(range2iono.lat.values)):
        # print(str(i) + '/' + str(len(range2iono.lat.values)))
        for j in range(len(range2iono.lon.values)):
            if ~np.isnan(incml.values[i, j]):
                # theta = float(np.radians(incml.values[i, j]))
                eledeg = float(90 - incml.values[i, j])
                ilat_ground, ilon_ground = range2iono.lat.values[i], range2iono.lon.values[j]
                x, y, z = aer2ecef(azimuthDeg, eledeg, range2iono.values[i, j], ilat_ground, ilon_ground, 0)
                                 #  float(hgtml.values[i, j])) # to consider hgt ... better without
                ilat, ilon, ialt = ecef2latlonhei(x, y, z)
                theta = float(np.radians(incml.values[i, j]))
                sin_thetaiono = earth_radius / (earth_radius + hiono) * np.sin(theta)
                if source=='code':
                    ionoij = get_vtec_from_tecxr(tecxr, acqtime, ilat, ilon)
                else:
                    ionoij = get_tecs(ilat, ilon, sat_alt_km, [acqtime], False)[0]
                ionoxr.values[i, j] = ionoij / np.sqrt(1 - sin_thetaiono ** 2) # with the last term, we get it to LOS (STEC)
    # now, convert TEC values into 'phase' - simplified here (?)
    f0 = 5.4050005e9
    # inc = avg_incidence_angle  # e.g. 39.1918 ... oh but... it actually should be the iono-squint-corrected angle. ignoring now
    # ionoxr = -4*np.pi*40.308193/speed_of_light/f0*ionoxr/np.cos(np.radians(incml))
    ionoxr = 4 * np.pi * 40.308193 / speed_of_light / f0 * ionoxr
    tecphase = ionoxr #get_tecphase(epoch)
    tecphase = interpolate_nans_bivariate(tecphase)
    tecphase = tecphase.interp_like(inc, method='linear', kwargs={"bounds_error": False, "fill_value": None})
    if np.max(np.isnan(tecphase.values)):
        tecphase = interpolate_nans_bivariate(tecphase)  # not the best (memory...) but needed
    return tecphase


# test frame: 144A_04689_111111
def make_all_frame_epochs(frame, source = 'code', epochslist = None):
    ''' use either 'code' or 'iri' as the source model for the correction
    Args:
        frame (str)
        source (str): either 'iri' or 'code'
        epochslist (list): e.g. ['20180930', '20181012'] - if given, only IPS for only those epochs are created, otherwise for all epochs
    '''
    framepubdir = os.path.join(os.environ['LiCSAR_public'], str(int(frame[:3])), frame)
    hgt = os.path.join(framepubdir, 'metadata', frame+'.geo.hgt.tif')
    hgt = load_tif2xr(hgt)
    mask = (hgt != 0) * (~np.isnan(hgt))
    if not epochslist:
        epochslist = os.listdir(os.path.join(framepubdir, 'epochs'))
    for epoch in epochslist:
        print(epoch)
        tif = os.path.join(framepubdir, 'epochs', epoch, epoch+'.geo.iono.'+source+'.tif')
        xrda = make_ionocorr_epoch(frame, epoch, source = source)
        xrda = xrda.where(mask)
        export_xr2tif(xrda, tif)

