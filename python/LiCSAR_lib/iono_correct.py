#!/usr/bin/env python3
import xarray as xr
from daz_iono import *
from lics_unwrap import *
from scipy.constants import speed_of_light
import numpy as np
from scipy.interpolate import griddata
from LiCSAR_misc import *
import framecare as fc

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
        orientation=frame[3]
        Efile = os.path.join(metadir, frame + '.geo.E.tif')
        Nfile = os.path.join(metadir, frame + '.geo.N.tif')
        E = load_tif2xr(Efile)
        N = load_tif2xr(Nfile)
        E = E.where(E != 0)
        N = N.where(N != 0)
        if orientation == "A":
            head_rad = np.arcsin(N / np.sin(np.radians(inc)))
            heading = np.degrees(head_rad)
        elif orientation == "D":
            head_rad = np.arcsin(- N / np.sin(np.radians(inc))) - np.pi
            heading = np.degrees(head_rad)   
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



def make_ionocorr_pair(frame, pair, sbovl=False,source = 'code', fixed_f2_height_km = 450, outif=None):
    """ This will generate ionospheric correction for given frame-pair.
    It would optionally output the result to a geotiff.
    
    Args:
        frame (str):    frame ID
        pair (str):     pair (e.g. '20180930_20181012')
        sbovl (bool):   if the sbovl true, ionospheric gradient calculated
        source (str):   source model for TEC values. Either 'iri' or 'code'.
        fixed_f2_height_km (int):  if None, it will estimate this using IRI
        outif (str):    if given, will export the iono phase screen to given geotiff
    Returns:
        xr.DataArray:   estimated ionospheric phase screen
    """
    # TODO - this below will mean we cannot use local ifgs. But can/must be (easily) improved!
    if not sbovl:
        ifg = load_ifg(frame, pair, unw = False, mag = False, prefer_unfiltered = False, stdout = False) # just to get mask etc.
        epochs = pair.split('_')
        tecphase1 = os.path.join(os.environ['LiCSAR_public'], str(int(frame[:3])),frame,'epochs',epochs[0],epochs[0]+'.geo.iono.'+source+'.tif')
        tecphase2 = os.path.join(os.environ['LiCSAR_public'], str(int(frame[:3])),frame,'epochs',epochs[1],epochs[1]+'.geo.iono.'+source+'.tif')
        if os.path.exists(tecphase1):
            tecphase1 = load_tif2xr(tecphase1)
        else:
            tecphase1 = make_ionocorr_epoch(frame, epochs[0], fixed_f2_height_km = fixed_f2_height_km, source = source)
        if os.path.exists(tecphase2):
            tecphase2 = load_tif2xr(tecphase2)
        else:
            tecphase2 = make_ionocorr_epoch(frame, epochs[1], fixed_f2_height_km = fixed_f2_height_km, source = source)
        # do their difference
        tecdiff = tecphase1 - tecphase2
        #    # tecdiff = interpolate_nans_pyinterp(tecdiff)
        #tecdiff = interpolate_nans_bivariate(tecdiff)
        #tecdiff = tecdiff.interp_like(ifg, method='linear', kwargs={"bounds_error": False, "fill_value": None})
        #tecdiff = interpolate_nans_bivariate(tecdiff) # not needed?
        #    if np.max(np.isnan(tecdiff.values)):
        #        tecdiff = interpolate_nans_bivariate(tecdiff)
        ifg['tecdiff'] = ifg['mask']*0.1
        ifg.tecdiff.values = tecdiff.values
        tecdiff = ifg.tecdiff.where(ifg.mask_extent == 1).rename('tecdiff_pha')
        if outif:
            export_xr2tif(tecdiff,outif)
        return tecdiff

    else:
        # print('sbovl iono correction')
        epochs = pair.split('_')
        tecsA1 = os.path.join(os.environ['LiCSAR_public'], str(int(frame[:3])),frame,'epochs',epochs[0],epochs[0]+'.geo.iono.'+source+'.sTECA.tif') #A,B=backward,forward lookings in azimuth; 1,2=prime,second epoch.
        tecsB1 = os.path.join(os.environ['LiCSAR_public'], str(int(frame[:3])),frame,'epochs',epochs[0],epochs[0]+'.geo.iono.'+source+'.sTECB.tif')
        tecsA2 = os.path.join(os.environ['LiCSAR_public'], str(int(frame[:3])),frame,'epochs',epochs[1],epochs[1]+'.geo.iono.'+source+'.sTECA.tif')
        tecsB2 = os.path.join(os.environ['LiCSAR_public'], str(int(frame[:3])),frame,'epochs',epochs[1],epochs[1]+'.geo.iono.'+source+'.sTECB.tif') 
        ##check the sTECs
        #epoch1
        if os.path.exists(tecsA1) and os.path.exists(tecsB1):
            tecsA1 = load_tif2xr(tecsA1)
            tecsB1 = load_tif2xr(tecsB1)
        else:
            tecsA1, tecsB1 = make_ionocorr_epoch(frame, epoch[0], source = source, fixed_f2_height_km = fixed_f2_height_km, alpha = alpha, sbovl = sbovl)
        #epoch2
        if os.path.exists(tecsA2) and os.path.exists(tecsB2):
            tecsA2 = load_tif2xr(tecsA2)
            tecsB2 = load_tif2xr(tecsB2)
        else:
            tecsA2, tecsB2 = make_ionocorr_epoch(frame, epoch[1], source = source, fixed_f2_height_km = fixed_f2_height_km, alpha = alpha, sbovl = sbovl)

        ###parameter for TEC gradient
        azpix=14000
        PRF = 486.486
        k = 40.308193 # m^3 / s^2
        f0 = 5.4050005e9
        c = speed_of_light
        ##scaling_tif
        metafolder = os.path.join(os.environ['LiCSAR_public'], str(int(frame[:3])), frame, 'metadata')
        # Check if the metadata folder exists
        if os.path.exists(metafolder) and os.path.isdir(metafolder):
            scaling_tif = None  # Initialize variable to track if a file is found
            
            for files in os.listdir(metafolder):  
                if files.endswith('.geo.sbovl_scaling.tif'):
                    scaling_tif = os.path.join(metafolder, files)
            # Check if no scaling file was found
            if scaling_tif is None:
                print("No .geo.sbovl_scaling.tif file found in metadata folder.")
        else:
            print(f"metadata is not exist in LiCSAR_public")  

        ##scaling2dfdc
        scaling_factor=load_tif2xr(scaling_tif)
        dfDC=azpix*PRF/(2*np.pi*scaling_factor)
        fH = f0 + dfDC*0.5
        fL = f0 - dfDC*0.5

        ##gradient method Lazecky et al. 2023,GRL #https://github.com/comet-licsar/daz/blob/main/lib/daz_iono.py#L561
        tecovl = (tecsA2 - tecsA1)/fH - (tecsB2 - tecsB1)/fL ##TODO ask the milan for confirmation: A2,B2 is secondary epoch and A1,B1 is primary epoch,  tecovl = (tec_A_master - selected_frame_esds['TECS_A'])/fH - (tec_B_master - selected_frame_esds['TECS_B'])/fL
        iono_grad = 2*PRF*k/c/dfDC * tecovl #unitless
        iono_grad_mm=iono_grad*azpix #mm
            
        if outif:
            export_xr2tif(iono_grad_mm,outif)
        return iono_grad_mm


def correct_iono_pair(frame, pair, ifgtype = 'diff_pha', dolocal = False, infile = None, source = 'code', fixed_f2_height_km = 450, outif=None):
    """ This will correct the ifg pair

    Args:
        frame (str):    frame ID
        pair (str):     pair (e.g. '20180930_20181012')
        ifgtype (str):  one of 'diff_pha', 'diff_unfiltered_pha', 'unw' , s(bovl) to correct (or any other extension but apart from 'unw' it will wrap the phase)
        dolocal (bool): if True, it will try to find the ifg in local GEOC folder
        infile (str):   if given, it will use this path to load the tif instead of loading from LiCSAR_public
        source (str):   source model for TEC values. Either 'iri' or 'code'.
        fixed_f2_height_km (int):  if None, it will estimate this using IRI
        outif (str):    if given, will export the iono phase screen to given geotiff
    Returns:
        xr.DataArray:  corrected ifg
    """
    sbovl=False
    
    if ifgtype == 'unw' or ifgtype == 'sbovldiff.adf.mm' or ifgtype == 'bovldiff.adf.mm':
        unw = True
        if ifgtype == 'sbovldiff.adf.mm' or ifgtype == 'bovldiff.adf.mm':
            sbovl=True
    else:
        unw = False
    if not infile:
        if ifgtype == 'diff_unfiltered_pha':
            prefer_unfiltered = True
        else:
            prefer_unfiltered = False
        ifgcube = load_ifg(frame, pair, unw = unw, dolocal = dolocal, mag = False, cliparea_geo = None,
                           prefer_unfiltered = prefer_unfiltered, stdout = False)
        if unw:
            ifg = ifgcube['unw']
        else:
            ifg = ifgcube['pha']
        ifg = ifg * ifgcube['mask']
        ifgcube = None
    else:
        ifg = load_tif2xr(infile)
    ifg = ifg.where(ifg != 0)
    tecdiff = make_ionocorr_pair(frame, pair, sbovl=sbovl, source=source, fixed_f2_height_km=fixed_f2_height_km, outif=None)
    ifg.values = ifg.values - tecdiff.values # AREA_OR_POINT might clash. Assuming same it may differ by 1/2 pixel (ok for ionosphere..)
    if not unw:
        ifg.values = wrap2phase(ifg.values)
    if outif:
        export_xr2tif(ifg, outif)
    return ifg


def make_ionocorr_epoch(frame, epoch, source = 'code', fixed_f2_height_km = 450, alpha = 0.85, return_phase = True, sbovl= False, outif=None):
    """
    Args:
        ...
        fixed_f2_height_km (int or None): CODE is valid for h=450. Still, if None, it will use IRI2016 to estimate the height (in mid point between scene centre and satellite)
        alpha (float): used only for 'CODE' - standard value is 0.85
        return_phase (bool): if not, it will return TEC, otherwise returns phase
        sbovl (bool): if True, it will calculate TEC values for BOI's different piercing points #TODO can be more precise for subswath overlap. 
        outif (str or None): output tif filename - if None, it will continue without saving

    Returns:
        xr.DataArray
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
    metafile = os.path.join(metadir,'metadata.txt')
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
    # centre_range_m=float(grep1line('centre_range_m',metafile).split('=')[1])
    try:
        centre_range_m = float(grep1line('centre_range_ok_m',metafile).split('=')[1]) # 2024: GAMMA had a bug wrongly informing on centre_range (may differ by 20 km or so!). fixed most of it
    except:
        centre_range_m=float(grep1line('centre_range_m',metafile).split('=')[1])
    #
    master=str(grep1line('master',metafile).split('=')[1])
    #
    acqtime = pd.to_datetime(str(epoch) + 'T' + center_time)
    # this is to get point between sat and scene centre
    theta = np.radians(avg_incidence_angle)
    wgs84 = nv.FrameE(name='WGS84')
    Pscene_center = wgs84.GeoPoint(latitude=scene_center_lat, longitude=scene_center_lon, degrees=True)
    burst_len = 7100*2.758277 #approx. satellite velocity on the ground 7100 [m/s] * burst_interval [s]
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
    path = nv.GeoPath(Pscene_center.to_nvector(), Psatg.to_nvector())   #https://pypi.org/project/nvector/
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
    ippg_lat, ippg_lon, ipp_alt = ecef2latlonhei(x, y, z)  ##IPP not used rest of the code? MN
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
    incml = inc.coarsen({'lat': mlfactorlat, 'lon': mlfactorlon}, boundary='trim').mean() ##downsampling prevents unnecesarry high-frequency noise.
    # get range towards iono single-layer in the path to the satellite, consider hgt
    # range2iono = (hiono - hgtml) / np.cos(np.radians(incml))
    # get range towards iono single-layer in the path to the satellite, do not consider hgt
    range2iono = hiono / np.cos(np.radians(incml))  ##IPP to ground in the slant range direction 
    earth_radius = 6378160  # m
    
    ##define the output xr
    if sbovl:
        ionoxrA = incml.copy(deep=True)
        ionoxrB = incml.copy(deep=True)
    else:
        ionoxr = incml.copy(deep=True)
    
    if source == 'code':
        tecxr=get_vtec_from_code(acqtime, lat=0, lon=0, return_fullxr = True)
        tecxr = alpha * tecxr
    print('getting TEC values sampled by {} km.'.format(str(round(ionosampling / 1000))))
    
    if not sbovl:  ##I could put the condition in the for loop but it's better to keep it clean.
        
        for i in range(len(range2iono.lat.values)):
            # print(str(i) + '/' + str(len(range2iono.lat.values)))
            for j in range(len(range2iono.lon.values)):
                if ~np.isnan(incml.values[i, j]):
                    # theta = float(np.radians(incml.values[i, j]))
                    eledeg = float(90 - incml.values[i, j])
                    #ground_scene
                    ilat_ground, ilon_ground = range2iono.lat.values[i], range2iono.lon.values[j]
                    
                    ##it directly starts from IPP scene
                    x, y, z = aer2ecef(azimuthDeg, eledeg, range2iono.values[i, j]*1000, ilat_ground, ilon_ground, 0)
                                    #  float(hgtml.values[i, j])) # to consider hgt ... better without
                    ilat, ilon, ialt = ecef2latlonhei(x, y, z)
                    
                    if source=='code':
                        ionoij = get_vtec_from_tecxr(tecxr, acqtime, ilat, ilon)
                    elif source=='iri':
                        ionoij = get_tecs(ilat, ilon, sat_alt_km, [acqtime], False)[0]  ##why sat_alt_km rather than ialt? MN
                        
                    else:
                        print(f'The source{source} model is not supported, try code or iri!')
                    
                    ##VTEC2STEC
                    theta = float(np.radians(incml.values[i, j]))
                    sin_thetaiono = earth_radius / (earth_radius + hiono) * np.sin(theta)
                    ionoxr.values[i, j] = ionoij / np.sqrt(1 - sin_thetaiono ** 2) # with the last term, we get it to LOS (STEC)
        
        if not return_phase:
            # if we want to return only TEC
            ionoxr = interpolate_nans_bivariate(ionoxr)
            ionoxr = ionoxr.interp_like(inc, method='linear', kwargs={"bounds_error": False, "fill_value": None})
            if np.max(np.isnan(ionoxr.values)):
                ionoxr = interpolate_nans_bivariate(ionoxr)  # not the best (memory...) but needed
            outputxr = ionoxr
        else:
            # now, convert TEC values into 'phase' - simplified here (?)
            f0 = 5.4050005e9
            # inc = avg_incidence_angle  # e.g. 39.1918 ... oh but... it actually should be the iono-squint-corrected angle. ignoring now
            # ionoxr = -4*np.pi*40.308193/speed_of_light/f0*ionoxr/np.cos(np.radians(incml))
            # we change the sign here as the effect on phase is opposite (ionospheric phase advance; less negative in case of stronger ionosphere)
            tecphase = 4 * np.pi * 40.308193 / speed_of_light / f0 * ionoxr
            # tecphase = ionoxr #get_tecphase(epoch)
            tecphase = interpolate_nans_bivariate(tecphase)
            tecphase = tecphase.interp_like(inc, method='linear', kwargs={"bounds_error": False, "fill_value": None})
            if np.max(np.isnan(tecphase.values)):
                tecphase = interpolate_nans_bivariate(tecphase)  # not the best (memory...) but needed
            outputxr = tecphase
        if outif:
            export_xr2tif(outputxr, outif)
        return outputxr           
                    
    elif sbovl:
        for i in range(len(range2iono.lat.values)):
            # print(str(i) + '/' + str(len(range2iono.lat.values)))
            for j in range(len(range2iono.lon.values)):
                if ~np.isnan(incml.values[i, j]):
                    # theta = float(np.radians(incml.values[i, j]))
                    eledeg = float(90 - incml.values[i, j])
                    #ground_scene
                    ilat_ground, ilon_ground = range2iono.lat.values[i], range2iono.lon.values[j]
                    ##satellite scene, we need to consider satellite scene again for BOI
                    x, y, z = aer2ecef(azimuthDeg, eledeg, slantRange, ilat_ground, ilon_ground, 0) #scene_alt) ### this is wrt ellipsoid! ##TODO rather than the center slant range near center and far range can be added! IDK it's effect.
                    satg_lat, satg_lon, sat_alt = ecef2latlonhei(x, y, z) ##satellite scene is changed becasue of the elevation degree change along the range direction
                    sat_alt_km = round(sat_alt / 1000)
                    Psatg = wgs84.GeoPoint(latitude=satg_lat, longitude=satg_lon, degrees=True)
                    #then get A', B' ##squint angle in azimuth direction lead that shift.
                    PsatgA, _azimuth = Psatg.displace(distance=burst_len/2, azimuth=heading-180, method='ellipsoid', degrees=True)
                    PsatgB, _azimuth = Psatg.displace(distance=burst_len/2, azimuth= heading, method='ellipsoid', degrees=True)
                    
                    ##IPP scene, 
                    x, y, z = aer2ecef(azimuthDeg, eledeg, range2iono.values[i, j]*1000, ilat_ground, ilon_ground, 0) #range should be in meter
                    ippg_lat, ippg_lon, ipp_alt = ecef2latlonhei(x, y, z)
                    Pippg = wgs84.GeoPoint(latitude=ippg_lat, longitude=ippg_lon, degrees=True)
                    #then get A', B'
                    PippAt, _azimuth = Pippg.displace(distance=burst_len, azimuth=heading-180, method='ellipsoid', degrees=True) #We extend burst_len to cover the intersection area. theorically burst_length/2 ideal
                    PippBt, _azimuth = Pippg.displace(distance=burst_len, azimuth= heading, method='ellipsoid', degrees=True)
        
                    ##intersection, to make sure IPP coordinates found
                    path_ipp = nv.GeoPath(PippAt, PippBt)
                    Pscene_center = wgs84.GeoPoint(latitude=ilat_ground, longitude=ilon_ground, degrees=True)
                    path_scene_satgA = nv.GeoPath(Pscene_center, PsatgA)
                    path_scene_satgB = nv.GeoPath(Pscene_center, PsatgB)
                    # these two points are the ones where we should get TEC
                    PippA = path_ipp.intersect(path_scene_satgA).to_geo_point()
                    PippB = path_ipp.intersect(path_scene_satgB).to_geo_point()
        
                    if source=='code':
                        # print('code selected')
                        ionoijA = get_vtec_from_tecxr(tecxr, acqtime, PippA.latitude_deg, PippA.longitude_deg)
                        ionoijB = get_vtec_from_tecxr(tecxr, acqtime, PippB.latitude_deg, PippB.longitude_deg)
                    else:
                        # print('iri selected')
                        ionoijA = get_tecs(PippA.latitude_deg, PippA.longitude_deg, sat_alt_km, [acqtime], False, source=source)[0]  ##TODO ask Milan why we set the sat_alt_km rather than ipp_alt when we apply IRI model? it select IPP itself?
                        ionoijB = get_tecs(PippB.latitude_deg, PippB.longitude_deg, sat_alt_km, [acqtime], False, source=source)[0]
        
                    #VTEC2STEC        
                    theta = float(np.radians(incml.values[i, j]))
                    sin_thetaiono = earth_radius / (earth_radius + hiono) * np.sin(theta)
                    ionoxrA.values[i, j] = ionoijA / np.sqrt(1 - sin_thetaiono ** 2)
                    ionoxrB.values[i, j] = ionoijB / np.sqrt(1 - sin_thetaiono ** 2)
        
        # Let's make it a bit more robust 
        ionoxr_dict = {'A': ionoxrA, 'B': ionoxrB}
        outputxr_dict = {}

        for key in ['A', 'B']:
            ionoxr_dict[key] = interpolate_nans_bivariate(ionoxr_dict[key])
            ionoxr_dict[key] = ionoxr_dict[key].interp_like(inc, method='linear', kwargs={"bounds_error": False, "fill_value": None})

            if np.max(np.isnan(ionoxr_dict[key].values)):
                ionoxr_dict[key] = interpolate_nans_bivariate(ionoxr_dict[key])  # Not the best, but needed

            outputxr_dict[key] = ionoxr_dict[key]

        # Return the processed values
        return outputxr_dict['A'], outputxr_dict['B']
                                    
                    
# test frame: 144A_04689_111111
def make_all_frame_epochs(frame, source = 'code', epochslist = None, fixed_f2_height_km = 450, alpha = 0.85, sbovl = False, startdate='20141001', enddate=None):
    ''' use either 'code' or 'iri' as the source model for the correction
    This function will generate ionosphere phase screens (LOS) [rad] per epoch.
    
    Args:
        frame (str)
        source (str): either 'iri' or 'code'
        epochslist (list): e.g. ['20180930', '20181012'] - if given, only IPS for only those epochs are created, otherwise for all epochs
        fixed_f2_height_km (num): for CODE only. CODE is valid for 450 km. if None, it will use IRI to estimate ionospheric height
        alpha (float): for CODE only
        sbovl (bool): if True, it will calculate TEC values for BOI's different piercing points
        startdate (str): Start date for filtering epochs, format "YYYYMMDD".
        enddate (str): End date for filtering epochs, format "YYYYMMDD".
    '''
    if enddate is None:
        enddate = pd.to_datetime("today").strftime("%Y%m%d")

    startdate = pd.to_datetime(startdate, format="%Y%m%d")
    enddate = pd.to_datetime(enddate, format="%Y%m%d")
    
    framepubdir = os.path.join(os.environ['LiCSAR_public'], str(int(frame[:3])), frame)
    hgtfile = os.path.join(framepubdir, 'metadata', frame+'.geo.hgt.tif')
    hgt = load_tif2xr(hgtfile)
    mask = (hgt != 0) * (~np.isnan(hgt))
    if not epochslist:
        epochslist = list(set(fc.get_epochs(frame) + fc.get_epochs_from_ifg_list_pubdir(frame)))
        epochslist.sort()
        #epochslist = os.listdir(os.path.join(framepubdir, 'epochs')) # careful, non-epoch folders would cause error!
    for epoch in epochslist:
        epoch_date = pd.to_datetime(epoch, format="%Y%m%d")
        
        if (epoch_date < startdate or epoch_date > enddate): ##just to scale between given time
            continue  # Skip epochs outside the selected date range if sbovl is True
 
        epochdir = os.path.join(framepubdir, 'epochs', epoch)
        if not os.path.exists(epochdir):
            os.mkdir(epochdir)
        
        if not sbovl:
            tif = os.path.join(epochdir, epoch+'.geo.iono.'+source+'.tif')
            if not os.path.exists(tif):
                print(epoch)
                xrda = make_ionocorr_epoch(frame, epoch, source = source, fixed_f2_height_km = fixed_f2_height_km, alpha = alpha, sbovl = sbovl)
                xrda = xrda.where(mask)
                export_xr2tif(xrda, tif) #, refto=hgtfile)
                os.system('chmod 777 '+tif)
                # but it still does not really fit - ok, because the xarray outputs here are gridline-registered while our ifgs are pixel registered...hmmm..
        elif sbovl:
            tif1= os.path.join(epochdir, epoch+'.geo.iono.'+source+'.sTECA.tif')
            tif2= os.path.join(epochdir, epoch+'.geo.iono.'+source+'.sTECB.tif')
            if not os.path.exists(tif1) or not os.path.exists(tif2):
                print(epoch)
                xrdaA, xrdaB = make_ionocorr_epoch(frame, epoch, source = source, fixed_f2_height_km = fixed_f2_height_km, alpha = alpha, sbovl = sbovl)
                #mask
                xrdaA = xrdaA.where(mask)
                xrdaB = xrdaB.where(mask)
                #export
                export_xr2tif(xrdaA, tif1) #, refto=hgtfile)
                export_xr2tif(xrdaB, tif2) #, refto=hgtfile)
                #permmssion
                os.system('chmod 777 '+tif1)
                os.system('chmod 777 '+tif2)
