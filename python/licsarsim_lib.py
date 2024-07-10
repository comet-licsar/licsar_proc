#!/usr/bin/env python3
# Simulation of SAR Intensity based on DEM
# M. Lazecky, 2024 (for Shailza Sharma, DEEPVOLC postdoc)

'''
Steps:
for one temporal epoch:
Inputs:
 - a DEM (geotiff)
 - heading, centre inc angle, centre inc angle
 - other params (corresponding to Sentinel-1 IW 'mid-swath')
Outputs:
 - geocoded simulated intensity

This will:
 - estimate satellite position (update SOVs) - primary contribution here
 - use gamma script to simulate S1 SAR intensity directly in geo coordinates (same resolution as input DEM)
    - this script can also output in rdc if needed

Naming convention of 
A) sample data:
e.g. 1000.054A.geo.mli.radcal.tif:
       1    2   3   4     5
 - 1 - volcano clip ID (LiCSVolc database)
 - 2 - relative orbit number + satellite pass (A for ascending, D for descending)
 - 3 - coordinates (geo: geographic system, rdc: radar coordinates)
 - 4 - type; mli = multilooked intensity
 - 5 - 'radiometrically calibrated'

B) sample outputs: (i.e. running the Example below on all DEM and related real data par files)
e.g. simsar.H-13.I39.1000.054A.geo.tif:
              1    2   3    4
 - 1 - applied heading (by default same as source mli.par file but can be changed for the simulation) [degrees from N]
 - 2 - applied incidence angle at the central pixel (see 1) [degrees, this angle is on the ground between the vertical and the satellite]
 - 3 - see A.1
 - 4 - see A.2
 

Versions:
0.0.1:
 - output is simulated SAR intensity geocoded in the same resolution as the input DEM
 
Example to generate for all available frames per given volcano:
volclip='23'
from licsarsim_lib import *
indem = volclip+'.dem'
for parfile in glob.glob(volclip+'.????.mli.par'):
    h,i,r = get_h_i_r_from_parfile(parfile)
    extraext = parfile[:-8]
    main_simsar(indem, h,i,r, extraext)
    
## parfile = '1000.054A.mli.par'
## h,i,r = get_h_i_r_from_parfile(parfile)
## extraext = parfile[:-8]  # extra text in output filenames
## indem = parfile.split('.')[0]+'.dem'
## main_simsar(indem, h,i,r, extraext)

# to preview:
from lics_vis import vis_tif; vis_tif('simsar.H-13.I39.1000.054A.geo.tif')
# to preview the orig (radiometrically calibrated) mli:
vis_tif('1000.163D.geo.mli.radcal.tif', to_amp_db = True)
# vis_tif('1000.163D.geo.mli.tif', to_amp_db = True)
# 
# NOTE: the input mlis should be calibrated first! 
'''

from orbit_lib import *
from daz_iono import *
import nvector as nv
# better to use py_gamma directly, but will just run direct commands
import os, glob
import subprocess as subp
from LiCSAR_misc import get_param_gamma
#from daz_lib_licsar import get_param_gamma   # daz was using framecare
import rioxarray
import xarray as xr
import numpy as np
import time

# INPUTS:
# DEM must be clipped to target area - you should use Geotiff (but gamma DEM should also work - just it might fail in some later step if par does not exist)
'''
# adding calibrated mlis:
pp=/gws/nopw/j04/nceo_geohazards_vol1/projects/LiCS/proc/current/subsets/volc
OUTDIR=/gws/nopw/j04/nceo_geohazards_vol1/projects/LiCS/proc/current/subsets/volc/for_simsar
cd $pp
for vv in `ls -d [1-9]*[0-9]`; do
 echo $vv
 cd $pp
 for fr in `ls $vv`; do
  cd $pp/$vv/$fr
  m=`ls SLC | head -n1`
  outname=$vv.$fr
  radcal_MLI SLC/$m/$m.slc.mli SLC/$m/$m.slc.mli.par - SLC/$m/$m.slc.mli.calibrated - 1 >/dev/null 2>/dev/null;
  mv SLC/$m/$m.slc.mli SLC/$m/$m.slc.mli.orig;
  cd SLC/$m; ln -s $m.slc.mli.calibrated $m.slc.mli; cd ../..;
  mv GEOC.MLI.30m/$m GEOC.MLI.30m/$m.uncalibrated;
  if [ ! -d geo ]; then ln -s geo.30m geo; fi;
  if [ ! -d GEOC.MLI ]; then ln -s GEOC.MLI.30m GEOC.MLI; fi;
  create_geoctiffs_to_pub.sh -M `pwd` $m >/dev/null 2>/dev/null; rm geo GEOC.MLI;
  cp GEOC.MLI.30m/$m/$m.geo.mli.tif $OUTDIR/$outname.geo.mli.radcal.tif
  # return it back
  rm SLC/$m/$m.slc.mli
  mv SLC/$m/$m.slc.mli.orig SLC/$m/$m.slc.mli
 done
done
'''


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


def main_simsar(indem, heading = -13, incidence_angle = 32, center_range_slc = 820000, extraext = ''):
    startime = time.time()
    simparams = extract_simparams(indem, heading, incidence_angle, center_range_slc)
    outfile = simulate_intensity(indem, simparams, extraext = extraext)
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


def simulate_intensity(indem = 'dem_crop.dem', simparams = None, extraext = ''):
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
    if extraext:
        strid = strid + '.' + extraext
    mlipar = 'simsar.'+strid+'.par'
    if not os.path.exists(mlipar):
        # prep some of the params:
        #  simparams = extract_simparams(dempar, simparams)
        generate_mlipar(mlipar, simparams)
    #
    # minimalistically to get only intensity:
    demseg = 'demseg'
    # this will be needed for parallelism.. so.. turning it on
    if extraext:
        demseg = demseg + '.' + extraext
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
        #print('note, simsar output is probably amplitude [dB], i.e. log10(sqrt(intensity))')
        print("from lics_vis import vis_tif; vis_tif('"+simsartif+"')")
    return simsartif




# 2024/06: for despeckling:
# 1, export MLIs to tif (in RDCs)
# 2, apply despeckling
# 3, resample to geocode

def rslc2mli(rslc, outtif = None):
    '''this will create intensity (i.e. amplitude squared) from the complex RSLC
    to convert the mli to dB, one may just do   10*np.log10(np.sqrt(mli))
    '''
    outmli = rslc[:-4]+'mli'
    if not outtif:
        outtif = outmli+'.tif'
    cmd = 'multi_look {0} {0}.par {1} {1}.par 1 1 >/dev/null'.format(rslc, outmli)
    os.system(cmd)
    wid = get_param_gamma('range_samples', outmli+'.par', floatt=False) # str
    cmd = 'data2tiff {0} {1} 2 {2} >/dev/null'.format(outmli, wid, outtif)
    os.system(cmd)
    # cleaning
    if os.path.exists(outmli):
        os.remove(outmli)
        #os.remove(outmli+'.par')
    return outtif


def geocode_tif(tif, lut = 'geo/20160901.lt_fine', 
                mlipar = 'RSLC/20220322/20220322.mli.par', 
                eqapar = 'geo/EQA.dem_par', diffpar = 'geo/20160901.diff_par', outtif = None):
    '''once despeckled and stored as tif in RDC, you should be able to geocode given a LUT
    Note, tif should be the despeckled tif, mlipar can be any previously existing mli.par file (must have same dimensions)'''
    a = rioxarray.open_rasterio(tif)
    tmpbin = tif+'.bin'
    if not outtif:
        outtif = tif[:-3]+'geo.tif'
    a.values.byteswap().astype(np.float32).tofile(tmpbin)
    width = len(a['x'])
    # first multilook
    rgfactor = get_param_gamma('range_looks', diffpar, floatt=False) # str
    azfactor = get_param_gamma('azimuth_looks', diffpar, floatt=False) # str
    cmd = 'multi_look_MLI {0} {1} {0}2 {0}2.par {2} {3} >/dev/null'.format(tmpbin, mlipar, rgfactor, azfactor)
    os.system(cmd)
    width_lut = get_param_gamma('width', eqapar, floatt=False) # str
    width_mli = get_param_gamma('range_samples',tmpbin+'2.par', floatt=False)
    cmd = 'geocode_back {0}2 {1} {2} {3} {4} - 1 0 >/dev/null'.format(tmpbin, width_mli, lut, outtif+'.bin', width_lut)
    print('geocoding')
    os.system(cmd)
    cmd = 'data2geotiff {0} {1} 2 {2} >/dev/null'.format(eqapar, outtif+'.bin', outtif)
    os.system(cmd)
    # cleaning
    if os.path.exists(outtif+'.bin'):
        os.remove(outtif+'.bin')
    if os.path.exists(tmpbin):
        os.remove(tmpbin)
    if os.path.exists(tmpbin+'2'):
        os.remove(tmpbin+'2')
    if os.path.exists(tmpbin+'2.par'):
        os.remove(tmpbin+'2.par')
    return outtif



def rslc2tif(rslc, outtif = None ):
    '''this will convert given SLC in binary format to complex-valued TIF
    '''
    if not outtif:
        outtif = rslc+'.tif'
    wid = get_param_gamma('range_samples', rslc+'.par', floatt=False) # str
    datatype = get_param_gamma('image_format', rslc + '.par', floatt=False)  # str
    if datatype == 'SCOMPLEX':
        gammatype = '3'
    elif datatype == 'FCOMPLEX':
        gammatype = '4'
    else:
        print('Something wrong with the data type - is this COMPLEX SLC?')
        return False
    cmd = 'data2tiff {0} {1} {2} {3} >/dev/null'.format(rslc, wid, gammatype, outtif)
    os.system(cmd)
    # this below is ugly but working
    a=rioxarray.open_rasterio(outtif)
    b=a.copy()
    b.values=a.real.astype(np.int16)
    c=a.copy()
    c.values=a.imag.astype(np.int16)
    c = c.assign_coords({'band': [2]})
    out = xr.concat([b,c], dim='band')
    os.system('rm '+outtif)
    out.rio.to_raster(outtif)
    return outtif

'''
a=rioxarray.open_rasterio('20160901.slc.tif')

b=a.copy()
b.values=a.real.astype(np.int16)

c=a.copy()
c.values=a.imag.astype(np.int16)
c = c.assign_coords({'band': [2]})
import xarray as xr

out = xr.concat([b,c], dim='band')
out.rio.to_raster('test2.tif')

'''


'''
inpath=/gws/nopw/j04/nceo_geohazards_vol1/projects/LiCS/proc/current/subsets/volc/20
for fr in `ls $inpath`; do
  echo $fr
  mkdir $fr
  cp -r $inpath/$fr/RSLC $fr/.  # i know.. but.. fast
  mkdir $fr/geo
  cp $inpath/$fr/geo*/EQA.dem_par $inpath/$fr/geo*/*.lt_fine $inpath/$fr/geo*/2*.diff_par $fr/geo/.
done

py:  this below should do whole job.. given the despeckle function is included J
import glob
fr='015A'
frames = ['015A', '095D', '117A', '168D']
for fr in frames:
    print(fr)
    rslcs = glob.glob(fr+'/RSLC/20*/*.rslc')
    lut = glob.glob(fr+'/geo/*lt_fine')[0]
    diffpar = glob.glob(fr+'/geo/*diff_par')[0]
    eqapar = fr+'/geo/EQA.dem_par'
    mlipar = None
    for rslc in rslcs:
        try:
            rdctif = rslc2mli(rslc)
            print(rdctif)
        except:
            print('erroneous')
        if not mlipar:
            mlipar = rslc[:-4]+'mli.par'
            os.system('cp '+mlipar+' '+fr+'/geo/.')
        # then step 2 despeckle...
        # rdctif = '/path/to/despeckled.tif'
        #geocode_tif(rdctif)



'''
