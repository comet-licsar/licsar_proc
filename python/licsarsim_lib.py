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

# INPUTS:
# must be clipped to target area - you should use Geotiff (but gamma DEM should also work - just it might fail in some later step if par does not exist)
#indem = 'dem_crop.dem'
heading = -13.5 # INPUT
incidence_angle = 32.2255 # INPUT
#range_inc = # no need?
center_range_slc = 820000 # INPUT
resolution_azi = 14.0 # INPUT
resolution_rg = 2.3 # INPUT
#sat_velocity = 7100 # from memory..
prf = 486.486 # Hz gives time per azi line (simplified?)

cenlon, cenlat, reslon, reslat, lenlon, lenlat = get_centre_lonlat_etc(demtif)
demres_approx_m = (reslon+reslat)/2 * 111111
mlazi = int(np.floor(demres_approx_m/resolution_azi)) # expecting same DEM resolution in both lon lat
mlrg = int(np.floor(demres_approx_m/resolution_rg))
spacing_rg = mlrg*resolution_rg
spacing_azi = mlazi*resolution_azi
# need to convert lenlon/lat to lines/samples using heading (and inc angle?) info
# but it might/should work with nan->0, i hope - need to work fast
nlines = lenlat
nsamples = lenlon
dist_az = spacing_azi * nlines # [m]
dist_rg = spacing_rg * nsamples # [m]

satpos_centre = get_sov_pos(cenlon, cenlat, heading, incidence_angle, center_range_slc)
# to get start/end ground coords (estimated)
wgs84 = nv.FrameE(name='WGS84')
centerpoint = wgs84.GeoPoint(latitude=cenlat, longitude=cenlon, z = 0, degrees=True)
startcoords, _azimuth = centerpoint.displace(distance=dist_az/2, azimuth=heading-180, method='ellipsoid', degrees=True)
endcoords, _azimuth = centerpoint.displace(distance=dist_az/2, azimuth=heading, method='ellipsoid', degrees=True)

satpos_start = get_sov_pos(startcoords.longitude_deg, startcoords.latitude_deg, heading, incidence_angle, center_range_slc)
satpos_end = get_sov_pos(endcoords.longitude_deg, endcoords.latitude_deg, heading, incidence_angle, center_range_slc)

#near and far range:
npoint, _azimuth = centerpoint.displace(distance=dist_rg/2, azimuth=heading-90, method='ellipsoid', degrees=True)
fpoint, _azimuth = centerpoint.displace(distance=dist_rg/2, azimuth=heading+90, method='ellipsoid', degrees=True)

def get_distance_ecef(ecef1, ecef2):
    # ecefs are tuples of (x,y,z)
    return np.sqrt( (ecef1[0]-ecef2[0])**2 + (ecef1[1]-ecef2[1])**2 + (ecef1[2]-ecef2[2])**2 )


nground_ecef = latlonhei2ecef(npoint.latitude_deg, npoint.longitude_deg, 0)
fground_ecef = latlonhei2ecef(fpoint.latitude_deg, fpoint.longitude_deg, 0)
nrange = get_distance_ecef(satpos_centre, nground_ecef)
frange = get_distance_ecef(satpos_centre, fground_ecef)

simparams = dict()
simparams['start_time'] = 1010.0
simparams['center_time'] = simparams['start_time'] + dist_az/resolution_azi/prf/2
simparams['end_time'] = simparams['start_time'] + dist_az/resolution_azi/prf
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
simparams['incidence_angle'] = incidence_angle # not used in processing
simparams['azimuth_angle'] = 90.0
simparams['radar_frequency'] = 5.4050005e+09
simparams['number_of_state_vectors'] = 3 # start, centre, end acquisition time
simparams['time_of_first_state_vector'] = simparams['start_time']
simparams['state_vector_interval'] = dist_az/resolution_azi/prf/2 # or 10
# then add:
simparams['state_vector_position_1'] = str(satpos_start[0])+'  '+str(satpos_start[1])+'  '+str(satpos_start[2]) 
simparams['state_vector_position_2'] = str(satpos_centre[0])+'  '+str(satpos_centre[1])+'  '+str(satpos_centre[2])
simparams['state_vector_position_3'] = str(satpos_end[0])+'  '+str(satpos_end[1])+'  '+str(satpos_end[2])
svi=simparams['state_vector_interval']
# keeping the velocity the same
simparams['state_vector_velocity_1'] = str((satpos_centre[0]-satpos_start[0])/svi)+'  '+str((satpos_centre[1]-satpos_start[1])/svi)+'  '+str((satpos_centre[2]-satpos_start[2])/svi)
simparams['state_vector_velocity_2'] = simparams['state_vector_velocity_1']
simparams['state_vector_velocity_3'] = simparams['state_vector_velocity_1']


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

'''
def runcmd(cmd, message):
    print(message)
    os.system(cmd + ' >/dev/null 2>/dev/null')
'''

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


#from LiCSAR_misc import grep1
from daz_lib_licsar import get_param_gamma
import rioxarray
import numpy as np
def get_centre_lonlat_etc(demtif):
    a = rioxarray.open_rasterio(demtif)
    lon = float(a['x'].mean())
    lat = float(a['y'].mean())
    reslon = np.abs(float(a.x[1]-a.x[0]))
    reslat = np.abs(float(a.y[1]-a.y[0]))
    lenlon = len(a.x)
    lenlat = len(a.y)
    return lon, lat, reslon, reslat, lenlon, lenlat


def extract_simparams(dempar, simparams):
    #get_param_gamma(param, parfile, floatt = True, pos = 0)
    
    center_latitude
    center_longitude
    indempar:
        corner_lat
        corner_lon
    grep1(arg,dempar)
    simparams[''] = get_param_gamma('corner_lat', dempar)
    return simparams


def simulate_intensity(indem = 'dem_crop.dem', simparams = simparams):
    # convert to gamma format (or use existing)
    if indem[:-3] == 'tif':
        dembin = indem[:-4]
        dodemtif = False
        demtif = indem
    else:
        dembin = indem
        dodemtif = True
    if os.path.exists(dembin+'_par'):
        dempar = dembin+'_par'
    else:
        dempar = dembin+'.par'
    
    geoid = '-' # or $DIFF_HOME/scripts/egm2008-5.dem or $DIFF_HOME/scripts/egm96.dem
    geoidpar = '-' # or:
    # geoidpar = geoid+'.dem_par'
    
    if not os.path.exists(dempar):
        cmd = ['dem_import', indemtif, dembin, dempar, '-', '-', geoid, geoidpar]
        runcmd(cmd, "Converting DEM to GAMMA format (geo)")
    
    # create needed files:
    mlipar = 'simsar.par'
    if not os.path.exists(mlipar):
        # prep some of the params:
        simparams = extract_simparams(dempar, simparams)
        generate_mlipar(mlipar, simparams)
    
    # minimalistically to get only intensity:
    demseg = 'demseg'
    demsegpar = demseg+'.par'
    lut = 'lut'
    lsmap = '-'
    incmap = '-'
    resmap = '-'
    lamap = '-'
    simsar = 'simsar.geo'   # output seems in dB (intensity->log10)
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
    runcmd(cmd, "Simulating DEM amplitude in RDC")
    
    # now to convert simsar to something normal, e.g.:
    gdtype = '2' # 2=FLOAT
    simsartif = simsar+'.tif'
    cmd = ['data2geotiff', demsegpar, simsar, gdtype, simsartif]
    cmdone = runcmd(cmd, "Exporting to "+simsartif)
    if dodemtif:
        demtif = demseg+'.tif'
        cmd = ['data2geotiff', demsegpar, demseg, gdtype, demtif]
        runcmd(cmd)
    #pixareamaptif = pixareamap+'.tif'
    #cmd = ['data2geotiff', demsegpar, pixareamap, gdtype, pixareamaptif]
    #runcmd(cmd, "Exporting to "+pixareamaptif)
    if cmdone:
        print('done. to preview, do (in python):')
        print("from lics_vis import vis_tif; vis_tif('simsar.geo.tif')")


# the below might not be needed:
# exit()
'''







outsig = 'sim_sigma0.rdc'
outgam = 'sim_gamma0.rdc'
# generating pix_sigma0 that should be used (but check also gamma0!)
#pixel_area(mlipar,demseg,lut,lsmap,inc,pixsigma,pixgamma,logfilename):
cmd = ["pixel_area", mlipar, dempar, dembin, lut, lsmap, incmap, outsig, outgam]
runcmd(cmd, "Simulating DEM amplitude in RDC")



def geocode_dem(masterslcdir,geodir,demdir,procdir,masterdate,outres = gc.outres, skip_fit = False):
    """ Geocodes the DEM to the master image geometry
    """
############################################################ Create radar coded DEM
    print("\nCreating radar coded DEM:")
    try:
        masterstr = masterdate.strftime('%Y%m%d')
    except:
        masterstr = masterdate
    logfilename = os.path.join(procdir,'log','gc_map_{0}.log'.format(masterstr))
    #create mli parameter file path, dem file path, dem parameter path
    mlipar = os.path.join(masterslcdir,masterstr+'.slc.mli.par')
    if not os.path.exists(mlipar):
        mlipar = os.path.join(masterslcdir,masterstr+'.rslc.mli.par')
    dem = os.path.join(demdir,'dem_crop.dem')
    dempar = os.path.join(demdir, 'dem_crop.dem_par')
    # Get dem res N and E (which are latitude and longitude resolutions?)
    with open(dempar, 'r') as f:
        for line in f:
            if 'post_lat' in line:
                demresN = float(line.split()[1])
            if 'post_lon' in line:
                demresE = float(line.split()[1])
    #calculate oversampling factors for lon/lat
    ovrfactN = str(-1.0*(demresN/outres))
    ovrfactE = str((demresE/outres))
    #create a dem segment file path? map segment?
    demseg = os.path.join(geodir,'EQA.dem')
    #look up table path?
    lut = os.path.join(geodir,masterstr+'.lt')
    #simulated SAR backscatter image
    simsar = os.path.join(geodir,masterstr+'.sim_sar')
    #zenith angle of surface normal vector n 
    u = os.path.join(geodir,'u')
    #orientation angle of n 
    v = os.path.join(geodir,'v')
    #local incidence angle
    inc = os.path.join(geodir,'inc')
    #projection angle
    psi = os.path.join(geodir,'psi')
    #pixel area normalization factor
    pix = os.path.join(geodir,'pix')
    #layover and shadow map
    lsmap = os.path.join(geodir,'ls_map')
############################################################ Calculate look up table
    print("Calculating look-up table...")
    #create a geocoded lookup table. Does this estimate demseg par?
    if not gc_map(mlipar,'-',dem,demseg,lut,ovrfactN,ovrfactE,simsar,
                                      u,v,inc,psi,pix,lsmap,'8','2','',
                                      logfilename):
        print("\nError:", file=sys.stderr)
        print("Something went wrong during the DEM lookup table creation.", file=sys.stderr)
        return 1 
############################################################ Simulate DEM amplitude
    mli2 = mlipar[:-4]
    if not skip_fit:
        #sigma 0 normalization area
        pixsigma = os.path.join(geodir,'pix_sigma0')
        #gamma 0 normalization area
        pixgamma = os.path.join(geodir,'pix_gamma0')
        logfilename = os.path.join(procdir,'log','pixel_area_{0}.log'.format(masterstr))
        print('Simulating DEM amplitude...')
        if not pixel_area(mlipar,demseg,lut,lsmap,inc,pixsigma,pixgamma,logfilename):
            # Error in amplitude simulation
            print("\nError:", file=sys.stderr)
            print("Something went wrong during the DEM amplitude simulation.", file=sys.stderr)
            return 2
    ############################################################ Create a Diff param file
    #    masterstr = masterdate.strftime('%Y%m%d')
        #the pixsigma file from the previous step
        mli1 = os.path.join(geodir,'pix_sigma0')
        #mli2 = mlipar[:-4] #The original master mli
        #Diff param file path
        diffpar = os.path.join(geodir,masterstr+'.diff_par')
        logfile = os.path.join(procdir,'log',
                               'create_diff_par_{0}.log'.format(masterstr))
        if not create_diff_par(mli2+'.par','-',diffpar,'1','0',logfile):
            print("\nError:", file=sys.stderr)
            print("Something went wrong during the DIFF parameter file creation", file=sys.stderr)
            return 3
    ############################################################ Estimate cross correlation
                                                    #offsets between sim and master original
        print('Estimating offsets...')
        #offset file
        offfile = os.path.join(geodir,masterstr+'.offs')
        #redundent variable?
        cofffile = os.path.join(geodir,'coffs')
        #patch cross-correlation
        ccpfile = os.path.join(geodir,masterstr+'.ccp')
        #text file containing offsets
        offsets = os.path.join(geodir,'offsets')
        logfile = os.path.join(procdir,'log',
                               'offset_pwrm_{0}.log'.format(masterstr))
        if not offset_pwrm(mli1,mli2,diffpar,offfile,ccpfile,'256','256',
                           offsets,'2','64','64','0.2',logfile):
            print("\nError:", file=sys.stderr)
            print("Something went wrong during the cross correlation"\
                                                                "offset estimation.", file=sys.stderr)
            return 3
    ############################################################ Fit offset function
        coffs = os.path.join(geodir,'coffs')
        logfile = os.path.join(procdir,'log',
                               'offset_fitm_{0}.log'.format(masterstr))
        print('Fitting offsets...')
        if not offset_fitm(offfile,ccpfile,diffpar,coffs,coffs+'ets','0.2','1',logfile):
            print("\nError:", file=sys.stderr)
            print("Something went wrong during the offset function fittin.", file=sys.stderr)
            return 3 
    ############################################################ Refine lookup table
        demwidth,demlength = get_dem_size(os.path.join(geodir,'EQA.dem_par')) # isn't this demseg?
        lutfine = os.path.join(geodir,masterstr+'.lt_fine')
        logfile = os.path.join(procdir,'log',
                               'gc_map_fine_{0}.log'.format(masterstr))
        if not gc_map_fine(lut,demwidth,diffpar,lutfine,'1',logfile):
            # Error in lookup table refining
            print("\nError:", file=sys.stderr)
            print("Something went wrong refining the lookup table.", file=sys.stderr)
            return 4
############################################################ create dem seg mli
                                                #by gecoding master mli to demseg
    else:
        print('skipping fitting of MLI to DEM to improve geocoding')
        lutfine = lut
        demwidth,demlength = get_dem_size(os.path.join(geodir,'EQA.dem_par'))
    
    width, length = get_mli_size(mli2+'.par')
    demsegmli = os.path.join(geodir,'EQA.{0}.slc.mli'.format(masterstr))
    logfile = os.path.join(procdir,'log',
                           'geocode_back_{0}.log'.format(masterstr))
    print('Geocoding...')
    if not geocode_back(mli2,width,lutfine,demsegmli,demwidth,
                                        demlength,'2','0',logfile):
        # Error in lookup table refining
        print("\nError:", file=sys.stderr)
        print("Something went wrong creating the DEM mli.", file=sys.stderr)
        return 5
############################################################ Create a height and 
                                    #intensisty ras of the mli2 data based on DEM
    logfile = os.path.join(procdir,'log',
                           'rashgt_{0}_hgt.log'.format(masterstr))
    if int(width)/1000 > 1:
        reducfac = str(int(width)/1000)
    else:
        reducfac = '1'
    if not rashgt(demseg,demsegmli,str(demwidth),'-','-','-',reducfac,reducfac,'500',logfile):
        print("\nError:", file=sys.stderr)
        print("Something went wrong creating the DEM mli sunraster file", file=sys.stderr)
        return 5
############################################################ Geocode (resample) DEM 
                                                            # -> mli2 coords
    hgtfile = os.path.join(geodir,masterstr+'.hgt')
    logfile = os.path.join(procdir,'log',
                           'geocode_{0}.log'.format(masterstr))
    if not geocode(lutfine,demseg,str(demwidth),hgtfile,str(width),str(length),'2','0',logfile):
        print("\nError:", file=sys.stderr)
        print("Something went wrong creating the master heightfile", file=sys.stderr)
        return 6
############################################################ Create a height and 
                                                            #intensisty ras of the
                                                            #mli2 data based on DEM
    if not rashgt(hgtfile,mli2,str(width),'-','-','-',reducfac,reducfac,'500',logfile):
        print("\nError:", file=sys.stderr)
        print("Something went wrong creating the DEM mli sunraster file", file=sys.stderr)
        return 6
    return 0
'''
