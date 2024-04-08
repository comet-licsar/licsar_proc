#!/usr/bin/env python3
# Simulation of SAR Intensity based on DEM
# M. Lazecky, 2024 (for Shailza Sharma, DEEPVOLC postdoc)

'''
Steps:
for one temporal epoch:
Inputs:
 - a DEM (geotiff) - will be simulated elsewhere
 - heading, centre inc angle, inc angle spread (e.g. +-5 deg) - target rg/azi sampling? (maybe not)
 - other params (corresponding to Sentinel-1 IW 'mid-swath')
Outputs:
 - geocoded simulated intensity

This will probably:
 - estimate satellite position (update SOVs)
 - convert DEM to radar coords
 - create geocoding tables for this DEM
 - simulate rdc intensity
 - geocode the simulated intensity
 

Versions:
0.0.1:
 - output is simulated SAR intensity geocoded in the same resolution as the input DEM
'''

from orbit_lib import *
from daz_iono import *
import nvector as nv
# better to use py_gamma directly, but will just run direct commands
import os
import subprocess as subp
#from LiCSAR_misc import grep1
from daz_lib_licsar import get_param_gamma
import rioxarray
import numpy as np
import time

# INPUTS:
# DEM must be clipped to target area - you should use Geotiff (but gamma DEM should also work - just it might fail in some later step if par does not exist)

def get_h_i_r_from_parfile(parfile):
    ''' Gets inputs to main_simsar function from given par file

    Args:
        parfile (str): path to the par file

    Returns:
        heading, inc_angle, center_range_distance
    '''
    heading = get_param_gamma('heading', parfile, floatt=True, pos=0)
    inc = get_param_gamma('incidence_angle', parfile, floatt=True, pos=0)
    crange = get_param_gamma('center_range_slc', parfile, floatt=True, pos=0)
    return heading, inc, crange

def main_simsar(indem, heading = -13, incidence_angle = 32, center_range_slc = 820000):
    startime = time.time()
    simparams = extract_simparams(indem, heading, incidence_angle, center_range_slc)
    outfile = simulate_intensity(indem, simparams)
    print('')
    timeitsec = time.time() - startime
    print('Finished in {0} seconds. Output file: {1}'.format(str(np.round(timeitsec, 2)), outfile))


def extract_simparams(indem, heading = -13, incidence_angle = 32, center_range_slc = 820000,
                      resolution_azi = 14.0, resolution_rg = 2.3, prf = 486.486):
    ''' Extracts simulation parameters based on DEM info and some basic info on satellite geometry

    Args:
        indem (str):       path to the input DEM (TIF prefered but GAMMA format also ok)
        heading (float):    e.g. -13 for ascending, -169 for descending pass
        incidence_angle (float):    average inc angle in the DEM scene
        center_range_slc (float):   slant range distance (centre of DEM scene towards the satellite, over 800 km for S1)
        resolution_azi (float)
        resolution_rg (float)
        prf (float):    pulse repetition frequency [Hz]

    Returns:
        dict: simulation parameters
    '''
    demtif, dembin, dempar = check_convert_dem(indem)
    # get_param_gamma(param, parfile, floatt = True, pos = 0)
    cenlon, cenlat, reslon, reslat, lenlon, lenlat = get_centre_lonlat_etc(demtif)
    demres_approx_m = (reslon + reslat) / 2 * 111111
    mlazi = int(np.floor(demres_approx_m / resolution_azi))  # expecting same DEM resolution in both lon lat
    mlrg = int(np.floor(demres_approx_m / resolution_rg))
    spacing_rg = mlrg * resolution_rg
    spacing_azi = mlazi * resolution_azi
    # need to convert lenlon/lat to lines/samples using heading (and inc angle?) info
    # but it might/should work with nan->0, i hope - need to work fast
    nlines = lenlat
    nsamples = lenlon
    dist_az = spacing_azi * nlines  # [m]
    dist_rg = spacing_rg * nsamples  # [m]
    #
    satpos_centre = get_sov_pos(cenlon, cenlat, heading, incidence_angle, center_range_slc)
    # to get start/end ground coords (estimated)
    wgs84 = nv.FrameE(name='WGS84')
    centerpoint = wgs84.GeoPoint(latitude=cenlat, longitude=cenlon, z=0, degrees=True)
    startcoords, _azimuth = centerpoint.displace(distance=dist_az / 2, azimuth=heading - 180, method='ellipsoid',
                                                 degrees=True)
    endcoords, _azimuth = centerpoint.displace(distance=dist_az / 2, azimuth=heading, method='ellipsoid', degrees=True)
    #
    satpos_start = get_sov_pos(startcoords.longitude_deg, startcoords.latitude_deg, heading, incidence_angle,
                               center_range_slc)
    satpos_end = get_sov_pos(endcoords.longitude_deg, endcoords.latitude_deg, heading, incidence_angle,
                             center_range_slc)
    #
    # near and far range:
    npoint, _azimuth = centerpoint.displace(distance=dist_rg / 2, azimuth=heading - 90, method='ellipsoid',
                                            degrees=True)
    fpoint, _azimuth = centerpoint.displace(distance=dist_rg / 2, azimuth=heading + 90, method='ellipsoid',
                                            degrees=True)
    #
    nground_ecef = latlonhei2ecef(npoint.latitude_deg, npoint.longitude_deg, 0)
    fground_ecef = latlonhei2ecef(fpoint.latitude_deg, fpoint.longitude_deg, 0)
    nrange = get_distance_ecef(satpos_centre, nground_ecef)
    frange = get_distance_ecef(satpos_centre, fground_ecef)
    #
    simparams = dict()
    simparams['start_time'] = 1010.0
    simparams['center_time'] = simparams['start_time'] + dist_az / resolution_azi / prf / 2
    simparams['end_time'] = simparams['start_time'] + dist_az / resolution_azi / prf
    simparams['range_samples'] = nsamples
    simparams['azimuth_lines'] = nlines
    simparams['range_looks'] = mlrg
    simparams['azimuth_looks'] = mlazi
    simparams['image_format'] = 'FLOAT'
    simparams['image_geometry'] = 'SLANT_RANGE'
    simparams['range_scale_factor'] = 1.0
    simparams['azimuth_scale_factor'] = 1.0
    simparams['azimuth_deskew'] = 'ON'
    simparams['center_latitude'] = cenlat
    simparams['center_longitude'] = cenlon
    simparams['heading'] = heading  # not used in processing
    simparams['range_pixel_spacing'] = spacing_rg
    simparams['azimuth_pixel_spacing'] = spacing_azi
    simparams['near_range_slc'] = nrange
    simparams['center_range_slc'] = center_range_slc
    simparams['far_range_slc'] = frange
    simparams['incidence_angle'] = incidence_angle  # not used in processing
    simparams['azimuth_angle'] = 90.0
    simparams['radar_frequency'] = 5.4050005e+09
    simparams['number_of_state_vectors'] = 3  # start, centre, end acquisition time
    simparams['time_of_first_state_vector'] = simparams['start_time']
    simparams['state_vector_interval'] = dist_az / resolution_azi / prf / 2  # or 10
    # then add:
    simparams['state_vector_position_1'] = str(satpos_start[0]) + '  ' + str(satpos_start[1]) + '  ' + str(
        satpos_start[2])
    simparams['state_vector_position_2'] = str(satpos_centre[0]) + '  ' + str(satpos_centre[1]) + '  ' + str(
        satpos_centre[2])
    simparams['state_vector_position_3'] = str(satpos_end[0]) + '  ' + str(satpos_end[1]) + '  ' + str(satpos_end[2])
    svi = simparams['state_vector_interval']
    # keeping the velocity the same
    simparams['state_vector_velocity_1'] = str((satpos_centre[0] - satpos_start[0]) / svi) + '  ' + str(
        (satpos_centre[1] - satpos_start[1]) / svi) + '  ' + str((satpos_centre[2] - satpos_start[2]) / svi)
    simparams['state_vector_velocity_2'] = simparams['state_vector_velocity_1']
    simparams['state_vector_velocity_3'] = simparams['state_vector_velocity_1']
    return simparams


def get_distance_ecef(ecef1, ecef2):
    # ecefs are tuples of (x,y,z)
    return np.sqrt( (ecef1[0]-ecef2[0])**2 + (ecef1[1]-ecef2[1])**2 + (ecef1[2]-ecef2[2])**2 )


def get_sov_pos(glon, glat, heading, incangle, srange):
    '''
    azimuthDeg = heading-90 #yes, azimuth is w.r.t. N (positive to E)
    elevationDeg = 90-inc_angle_avg
    slantRange = range_avg
    '''
    #wgs84 = nv.FrameE(name='WGS84')
    #point = wgs84.GeoPoint(latitude=glat, longitude=glon, z = 0, degrees=True)
    sat_ecef = aer2ecef(heading-90, 90-incangle, srange, glat, glon, 0) # tuple
    return sat_ecef #... sat_ecef based


def get_resolution_tif(ifg, in_m=True):
    """Gets resolution of the xr.dataset (or dataarray), either in metres or degrees
    """
    resdeg = (np.abs(ifg.lat[1]-ifg.lat[0])+np.abs(ifg.lon[1]-ifg.lon[0]))/2
    if in_m:
        latres = 111.32 * np.cos(np.radians(ifg.lat.mean())) * 1000 # in m
        return float(latres * resdeg)
    else:
        return float(resdeg)


def runcmd(cmd, message = '',logdir = 'logs'):
    """ cmd can be either full string command or split to list (as for subp)
    """
    if message:
        print(message)
    if not os.path.exists(logdir):
        os.mkdir(logdir)
    if not isinstance(cmd, list):
        cmd = cmd.split()
    cmd0 = cmd[0]
    logfile = os.path.join(logdir, cmd0+'.log')
    errfile = os.path.join(logdir, cmd0+'.err')
    with open(logfile,'w') as f:
        with open(errfile,'w') as erf:
            try:
                rc = subp.check_call(cmd,stdout=f,stderr=erf)
            except:
                rc = 99
    if rc != 0:
        print('Something went wrong, see {0}'.format(logfile))
        return False
    else:
        return True


def generate_mlipar(parfile, simparams):
    if os.path.exists(parfile):
        os.remove(parfile)
    f = open(parfile, 'w')
    for param in simparams:
        val = simparams[param]
        f.write(param+':   '+str(val)+' \n')
    f.close()
    return


def get_centre_lonlat_etc(demtif):
    a = rioxarray.open_rasterio(demtif)
    lon = float(a['x'].mean())
    lat = float(a['y'].mean())
    reslon = np.abs(float(a.x[1]-a.x[0]))
    reslat = np.abs(float(a.y[1]-a.y[0]))
    lenlon = len(a.x)
    lenlat = len(a.y)
    return lon, lat, reslon, reslat, lenlon, lenlat


def check_convert_dem(indem, fix_geoid = False):
    if not fix_geoid:
        geoid = '-'
        geoidpar = '-'
    else:
        print('not implemented yet, see help of dem_import - need to add  $DIFF_HOME/scripts/egm2008-5.dem or $DIFF_HOME/scripts/egm96.dem')
        geoid = '-'
        geoidpar = '-'
    if indem[-3:] == 'tif':
        dembin = indem[:-4]
        dempar = dembin + '_par'
        if not os.path.exists(dempar):
            cmd = ['dem_import', indem, dembin, dempar, '-', '-', geoid, geoidpar]
            runcmd(cmd, "Converting DEM to GAMMA format (geo)")
        demtif = indem
    else:
        dembin = indem
        dempar = dembin + '_par'
        demtif = dembin+'.tif'
        if (not os.path.exists(demtif)) and (os.path.exists(dempar)):
            cmd = ['data2geotiff', dempar, dembin, '2', demtif]
            runcmd(cmd)
    return demtif, dembin, dempar


def simulate_intensity(indem = 'dem_crop.dem', simparams = None):
    ''' function to use simparams with the DEM to generate the simsar output

    Args:
        indem (str): path to the input DEM (should be tif but would work if in gamma format)
        simparams (dict): output of extract_simparams()

    Returns:
        str (path to the generated sim sar tiff)
    '''
    demtif, dembin, dempar = check_convert_dem(indem)
    #
    strid = 'H'+str(int(np.round(simparams['heading'])))+'.I'+str(int(np.round(simparams['incidence_angle'])))
    mlipar = 'simsar.'+strid+'.par'
    if not os.path.exists(mlipar):
        # prep some of the params:
        #  simparams = extract_simparams(dempar, simparams)
        generate_mlipar(mlipar, simparams)
    #
    # minimalistically to get only intensity:
    demseg = 'demseg'
    demsegpar = demseg+'.par'
    lut = 'lut.'+strid
    lsmap = '-'
    incmap = '-'
    resmap = '-'
    lamap = '-'
    simsar = 'simsar.'+strid+'.geo'   # output seems in dB (intensity->log10)
    #pixareamap = 'pixelarea'   # the normalisation gets too far from the SAR intensity!
    pixareamap = '-'
    '''
    # now get the oversampling right so the transformed DEM will have same dimensions as requested in mli.par
    # Get dem res N and E (which are latitude and longitude resolutions)
    outres = # get this from the mli par? but this would be azi/rg!
    with open(dempar, 'r') as f:
        for line in f:
            if 'post_lat' in line:
                demresN = float(line.split()[1])
            if 'post_lon' in line:
                demresE = float(line.split()[1])
    #calculate oversampling factors for lon/lat
    ovrfactN = str(-1.0*(demresN/outres))
    ovrfactE = str((demresE/outres))
    '''
    cmd = ['gc_map2', mlipar, dempar, dembin, demsegpar, demseg, lut, '-', '-', lsmap,'-', incmap, resmap, lamap, simsar, '-', '-', '-', pixareamap]
    runcmd(cmd, "Simulating DEM amplitude using gc_map2")
    #
    # now to convert simsar to something normal, e.g.:
    gdtype = '2' # 2=FLOAT
    simsartif = simsar+'.tif'
    cmd = ['data2geotiff', demsegpar, simsar, gdtype, simsartif]
    cmdone = runcmd(cmd, "Exporting to "+simsartif)
    #pixareamaptif = pixareamap+'.tif'
    #cmd = ['data2geotiff', demsegpar, pixareamap, gdtype, pixareamaptif]
    #runcmd(cmd, "Exporting to "+pixareamaptif)
    if cmdone:
        print('done. to preview, do (in python):')
        print("from lics_vis import vis_tif; vis_tif('simsar.geo.tif')")
    return simsartif
