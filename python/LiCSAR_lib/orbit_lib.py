
################################################################################
# Imports
################################################################################

import requests
import re
import datetime as dt
import os, shutil
import logging
import pandas as pd
from configparser import ConfigParser
import global_config as gc
from datetime import datetime
import xml.etree.ElementTree as ET
import pyproj
import numpy as np
import glob 

from sentinelsat.sentinel import SentinelAPI
try:
    from eof.download import download_eofs
except:
    print('sentineleof library is not installed (useful, see: https://github.com/scottstanie/sentineleof ) - no worries, we have workaround')
    #print('warning, you do not have sentineleof python library installed. expect problems with orbits (but maybe fixed now..)')


try:
    import nvector as nv
except:
    print('warning, nvector is not installed - some advanced functions will not work (but perhaps you dont need them)')

#### 2023 functions
#try:
#    from eof.parsing import parse_orbit
#except:
#    print('also, without sentineleof lib, we will not parse the orbits now')


# adapted from SentinelEOF parser:
def parse_utc_string(timestring):
    return datetime.strptime(timestring, "UTC=%Y-%m-%dT%H:%M:%S.%f")


# adapted from SentinelEOF parser:
def _convert_osv_field(osv, field, converter=float):
    # osv is a xml.etree.ElementTree.Element
    field_str = osv.find(field).text
    return converter(field_str)


# adapted from (useful) SentinelEOF parser:
def load_eof(
    eof_filename,
    min_time=datetime(1900, 1, 1),
    max_time=datetime(2100, 1, 1)):
    """ Would load the orbit into xr.Dataset
    Args:
        eof_filename: path to the orbits file
    Returns:
        xr.Dataset: time in UTC, coords in ECEF [m or m/s]
    """
    #
    tree = ET.parse(eof_filename)
    root = tree.getroot()
    all_osvs = []
    idxs_in_range = []
    for idx, osv in enumerate(root.findall("./Data_Block/List_of_OSVs/OSV")):
        all_osvs.append(osv)
        utc_dt = _convert_osv_field(osv, "UTC", parse_utc_string)
        if utc_dt >= min_time and utc_dt <= max_time:
            idxs_in_range.append(idx)
    #
    if not idxs_in_range:
        return None
    #
    idxs_in_range.sort()
    osvs_in_range = []
    for idx in idxs_in_range:
        cur_osv = all_osvs[idx]
        utc_dt = _convert_osv_field(cur_osv, "UTC", parse_utc_string)
        # utc_secs = secs_since_midnight(utc_dt)
        # cur_line = [utc_secs]
        cur_line = [utc_dt]
        for field in ("X", "Y", "Z", "VX", "VY", "VZ"):
            # Note: the 'unit' would be elem.attrib['unit']
            cur_line.append(_convert_osv_field(cur_osv, field, float))
        osvs_in_range.append(cur_line)
    #
    cols=['time','x','y','z','vx','vy','vz']
    osvs_in_range = pd.DataFrame(osvs_in_range, columns=cols)
    # convert to xarray (easy for later interpolations):
    orbxr = osvs_in_range.set_index('time').to_xarray()
    return orbxr


def get_coords_in_time(orbxr, timesample, method='cubic', return_as_nv = False):
    """ gets interpolated coordinates from the orbit datacube for given time sample (dt.datetime)
    Note, we need to implement hermite interpolation, as in ISCE2!!

    Args:
        orbxr (xr.Dataset):  e.g. using load_eof
        timesample (dt.datetime)
        method (str):  interpolation method (ML: note, linear vs cubic diff would reach tens of metres!)
    Returns:
        xr.Dataset
    """
    coords = orbxr.interp(time=timesample, method=method)
    if return_as_nv:
        lonlath1 = ecef2lonlathei(float(coords['x']), float(coords['y']),
                                  float(coords['z']))  # refElp.xyz_to_llh(vec1.getPosition())
        # import nvector as nv
        wgs84 = nv.FrameE(name='WGS84')
        point = wgs84.GeoPoint(latitude=lonlath1[1], longitude=lonlath1[0], z = lonlath1[2], degrees=True)
        return point
    else:
        return coords


def getoldorbpath(orbfiles):
    """helper function to get old orbit files corresponding to the given list"""
    oldorbs = []
    for orbfile in orbfiles:
        ff = os.path.basename(orbfile)
        #orbfile ='/gws/nopw/j04/nceo_geohazards_vol1/orbits_2021/S1A/POEORB/S1A_OPER_AUX_POEORB_OPOD_20210309T002908_V20180831T225942_20180902T005942.EOF'
        oldpath = '/gws/nopw/j04/nceo_geohazards_vol1/orbits.old/S1'+ff[2]+'/POEORB'
        oldorbcands = glob.glob(os.path.join(oldpath, 'S1'+ff[2]+'_OPER_AUX_POEORB_OPOD_*_V'+ff.split('V')[1]))
        oldorb = None
        for oldorbcand in oldorbcands:
            if not os.path.basename(oldorbcand) == ff:
                oldorb = oldorbcand
                break
        oldorbs.append(oldorb)
    return oldorbs


# from daz/daz_iono:
def ecef2lonlathei(x, y, z):
    transformer = pyproj.Transformer.from_crs(
        {"proj":'geocent', "ellps":'WGS84', "datum":'WGS84'},
        {"proj":'latlong', "ellps":'WGS84', "datum":'WGS84'},
        )
    lon, lat, alt = transformer.transform(x,y,z,radians=False)
    return lon, lat, alt


# should be right:
def get_azi_diff_from_two_orbits(orbfile1, orbfile2, timesample):
    """ Should get azi offset [m] from two different orbits.

    Args:
        orbfile1 (str): path to orbfile 1 (e.g. old orbit EOF file)
        orbfile2 (str): path to e.g. new orbit EOF
        timesample (dt.datetime)

    Returns:
        (float): azi offset [m]
    """
    oldorbxr = load_eof(orbfile1)
    neworbxr = load_eof(orbfile2)
    pointold = get_coords_in_time(oldorbxr, timesample, method='cubic', return_as_nv = True)
    pointoldPre = get_coords_in_time(oldorbxr, timesample-dt.timedelta(seconds=1), method='cubic', return_as_nv = True)
    pointnew = get_coords_in_time(neworbxr, timesample, method='cubic', return_as_nv = True)
    # now get azimuth direction shift
    heading = getHeading(neworbxr, time=timesample) # in degrees
    # get diff in azi
    pathA = nv.GeoPath(pointoldPre, pointold)
    pointC = pathA.closest_point_on_great_circle(pointnew)
    azidiff, _azi1, _azi2 = pointold.distance_and_azimuth(pointC)  # but check sign!!!
    
    # fix the sign
    azidiff = np.sign(_azi1)*np.sign(heading)*azidiff
    
    # get diff in rg
    #rgdiff = np.nan
    return azidiff #, rgdiff


def getHeading(orbxr, time=None, spacing=0.5):
    """
    Compute heading at given time.
    If time is not provided, mid point of orbit is used.
    Args:
        orbxr (xr.DataArray): datacube from state orbit vectors
        time (dt.datetime or None)
        spacing (float): timedelta to use for heading calculation, in seconds
    Returns:
        float: heading in degrees
    """

    if time is None:
        delta = orbxr.time.values.max() - orbxr.time.values.min()
        aztime = orbxr.time.values.min() + delta/2
        aztime = pd.to_datetime(str(aztime))
    else:
        aztime = time

    t1 = aztime - dt.timedelta(seconds=spacing)
    t2 = aztime + dt.timedelta(seconds=spacing)

    vec1 = get_coords_in_time(orbxr, t1)
    vec2 = get_coords_in_time(orbxr, t2)

    lonlath1 = ecef2lonlathei(float(vec1['x']), float(vec1['y']), float(vec1['z'])) #refElp.xyz_to_llh(vec1.getPosition())
    lonlath2 = ecef2lonlathei(float(vec2['x']), float(vec2['y']), float(vec2['z']))

    # see https://www.ffi.no/en/research/n-vector/#example_1
    # import nvector as nv
    wgs84 = nv.FrameE(name='WGS84')
    pointA = wgs84.GeoPoint(latitude=lonlath1[1], longitude=lonlath1[0], degrees=True)
    pointB = wgs84.GeoPoint(latitude=lonlath2[1], longitude=lonlath2[0], degrees=True)
    p_AB_N = pointA.delta_to(pointB)
    heading = p_AB_N.azimuth_deg
    return heading

'''
NOTES towards orb update:
timesample=dt.datetime(2016,12,14,10,10,45)
zipFile='S1B_IW_SLC__1SSV_20161214T100950_20161214T101017_003390_005CA6_6694'
orbitfiles = findValidOrbFile(baseDir,sat,startTime,endTime)

orbdir = os.environ['ORB_DIR']
orb='S1B_OPER_AUX_RESORB_OPOD_20161214T124022_V20161214T083137_20161214T114907'
eof_filename = getValidOrbFile(orbdir+'/S1B',orb)
neworbxr = load_eof(eof_filename)


oldorbxr = 

pointold = get_coords_in_time(oldorbxr, timesample, method='cubic', return_as_nv = True)
pointnew = get_coords_in_time(neworbxr, timesample, method='cubic', return_as_nv = True)

# now get azimuth direction shift
azimuthdir = getHeading(neworbxr, time=timesample) # in degrees

'''



'''
# adapted from isce2:
# https://github.com/isce-framework/isce2/blob/0dbb1679b61b4ac385c537269b91208e57024672/components/isceobj/Orbit/Orbit.py
def getHeading(orbxr, time=None, spacing=0.5):
    """
    Compute heading at given time.
    If time is not provided, mid point of orbit is used.
    Args:
        orbxr (xr.DataArray): datacube from state orbit vectors
        time (dt.datetime or None)
        spacing (float): timedelta to use for heading calculation, in seconds
    Returns:
        float: heading in degrees
    """

    if time is None:
        delta = orbxr.time.max() - orbxr.time.min()
        aztime = orbxr.time.min() + dt.timedelta(seconds = 0.5 * delta.total_seconds())
    else:
        aztime = time

    t1 = aztime - dt.timedelta(seconds=spacing)
    t2 = aztime + dt.timedelta(seconds=spacing)

    vec1 = get_coords_in_time(orbxr, t1)
    vec2 = get_coords_in_time(orbxr, t2)

    lonlath1 = ecef2lonlathei(float(vec1['x']), float(vec1['y']), float(vec1['z'])) #refElp.xyz_to_llh(vec1.getPosition())
    lonlath2 = ecef2lonlathei(float(vec2['x']), float(vec2['y']), float(vec2['z']))

    import nvector as nv
    wgs84 = nv.FrameE(name='WGS84')
    pointA = wgs84.GeoPoint(latitude=lonlath1[1], longitude=lonlath1[0], degrees=True)
    pointB = wgs84.GeoPoint(latitude=lonlath2[1], longitude=lonlath2[0], degrees=True)
    path = nv.GeoPath(pointA.to_nvector(), pointB.to_nvector())

    # TODO: can we do wgs84.point? and how to return the heading in degrees? AND.. once i am here, adapt also for diff between 2 orbit files in azimuth!
    #Heading
    return path.something
    #hdg = refElp.geo_hdg(lonlath1, lonlath2)
#
#    return np.degrees(hdg)
'''

'''
# from Reza B., to read state orbit vectors from S1 xmls (and orbit files).
# ML: to get diff between orbit files, I would:
# - load the SOVs to xarray
# - interpolate (cubic?) to get SOV for given time
# - do diff in y (should be azimuth?), x (range?) - or are the x,y,z ECEF coordinates instead of the satellite coords?


import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
from datetime import datetime
from scipy import interpolate
from astropy.time import Time
import framecare as fc
from mpl_toolkits.mplot3d import Axes3D

# Working directory
xml_dir = r'/home/users/eerbs/ion/20181011/'

# Path to the Sentinel-1 XML file
iw1_xml = xml_dir + 's1a-iw1-slc-hh-20181011t080210-20181011t080235-024085-02a1f4-001.xml'
iw2_xml = xml_dir + 's1a-iw2-slc-hh-20181011t080208-20181011t080233-024085-02a1f4-002.xml'
iw3_xml = xml_dir + 's1a-iw3-slc-hh-20181011t080209-20181011t080234-024085-02a1f4-003.xml'

# Parse the XML file and get roots
tree1 = ET.parse(iw1_xml)
tree2 = ET.parse(iw2_xml)
tree3 = ET.parse(iw3_xml)

root1 = tree1.getroot()
root2 = tree2.getroot()
root3 = tree3.getroot()

# State vectors extraction (satellite position and acquisition times) 
state_vectors = root1.findall('.//orbit')

pos_array = []
stateTime_array = []

for pos in state_vectors:
    x = float(pos.find('position/x').text)
    y = float(pos.find('position/y').text)
    z = float(pos.find('position/z').text)
    ts = pos.find('time').text

    stateTime_array.append(ts)
    pos_array.append([x,y,z])
'''


################################################################################
# Setup logger
################################################################################
logger = logging.getLogger(__name__)

################################################################################
# Exceptions
################################################################################
class badFileName(Exception):
    """Two options here - your filename doesn't contain the info we need or the pattern is broken"""
    def __init__(self,filename,timePat):
        self.filename = filename
        self.timePat = timePat
        
    def __str__(self):
        return "Could not parse {0} with pattern {1}".format(self.filename,self.timePat)

class invalidProductType(Exception):
    def __init__(self,prodType):
        self.prodType = prodType

    def __str__(self):
        return "Unregonised orbit product type {0}, use either POEORB or RESORB".format(self.prodType)

class invalidSatType(Exception):
    def __init__(self,satType):
        self.satType = satType

    def __str__(self):
        return "Unregonised satalite type {0}, use either S1A or S1B".format(self.satType)

class invalidFileUrl(Exception):
    """This means your url is probably wrong... where did you get it?"""
    def __init__(self,url):
        self.prodType = url

    def __str__(self):
        return "Could not parse an orbit file from url {0}".format(self.url)

class urlDownloadFailed(Exception):
    """Url download failed. Could just be internet connection issue"""
    def __init__(self,url):
        self.prodType = url

    def __str__(self):
        return "Could not download {0}".format(self.url)

class orbUrlNotFnd(Exception):
    """ Orbit file could not be found on the remote. Not always an issue"""
    def __init__(self,prodType,startTime,endTime):
        self.prodType = prodType
        self.startTime = startTime
        self.endTime = endTime
    def __str__(self):
        return "Could not find orbit url for product type {0} between times {1} and {2}".format(
                self.prodType,self.startTime,self.endTime)

################################################################################
# Strip satalite
################################################################################
def strpSat(filename):
    """ Strips out the satalite from the file name"""
    satPat = '.*(S1[AB])_.*'

    mtch = re.search(satPat,filename)

    if mtch:
        sat = mtch.groups()[0]
    else:
        raise badFileName(filename,satPat)

    return sat
################################################################################
# Strip file name to date times
################################################################################
def strpStrtEndTimes(filename):
    """ Strips out the measurement times from the file name"""
    timePat = '.*(\d{8}T\d+)_(\d{8}T\d+).*'
    strpPat = '%Y%m%dT%H%M%S'

    mtch = re.search(timePat,filename)

    if mtch:
        startTime = dt.datetime.strptime(mtch.groups()[0],strpPat)
        endTime = dt.datetime.strptime(mtch.groups()[1],strpPat)
    else:
        raise badFileName(filename,timePat)

    return (startTime,endTime)


def get_orbit_filenames_for_datetime(ddatetime, producttype='POEORB', s1ab = None):
    """
    gets orbits existing

    Args:
        ddatetime: dt.datetime or dt.datetime.date
        producttype: 'POEORB' or 'RESORB'
        s1ab (str): 'S1A' or 'S1B' - if None, it would return both
    Returns:
        list of filenames (full paths)
    """
    try:
        ddate = ddatetime.date()
    except:
        ddate = ddatetime
    listfiles = downloadOrbits_CopCloud(ddate-dt.timedelta(days=1), ddate+dt.timedelta(days=1), producttype)
    if s1ab:
        listf2 = []
        for f in listfiles:
            if os.path.basename(f)[:3] == s1ab:
                listf2.append(f)
        listfiles = listf2
    listf2 = []
    for f in listfiles:
        ff = os.path.basename(f)
        datein = ff.split('_')[6][1:]
        dateout = ff.split('_')[7].split('.')[0]
        if pd.Timestamp(datein)<ddatetime and pd.Timestamp(dateout)>ddatetime:
            listf2.append(f)
    listfiles = listf2
    return listfiles


def downloadOrbits_CopCloud(startdate, enddate, producttype):
    scihub = SentinelAPI('gnssguest', 'gnssguest','https://scihub.copernicus.eu/gnss')
    # for ONLY orbit files reprocessed in 2021
    result = scihub.query(platformname = 'Sentinel-1', producttype='AUX_'+producttype, date = (startdate, enddate), ingestionDate='[2021-01-01T00:00:00.000Z TO NOW]')
    # for 'any' orbit files
    #result = scihub.query(platformname = 'Sentinel-1', producttype='AUX_'+producttype, date = (startdate, enddate))    
    result = scihub.to_dataframe(result)
    existing = []
    for id, row in result.iterrows():
        outfile= os.path.join(os.environ['ORB_DIR'],'S1'+row['platformnumber'],producttype,row['filename'])
        if not os.path.exists(outfile):
            lockfile = outfile+'.lock'
            if os.path.exists(lockfile):
                print('orbit file locked - perhaps another orb download process running')
                continue
            else:
                f = open(lockfile, 'wb').close()
            #download it here....
            downurl = row.link.replace('https://','https://gnssguest:gnssguest@')
            try:
                print('downloading orbit file '+row.filename)
                r = requests.get(downurl, allow_redirects=True)
                f = open(outfile, 'wb')
                fsize = f.write(r.content)
                f.close()
            except:
                print('error downloading orbit file '+row.filename)
                print('trying from ASF - using wget')
                try:
                    parser = ConfigParser()
                    parser.read(gc.configfile)
                    asfuser = parser.get('asf', 'asfuser')
                    asfpass = parser.get('asf', 'asfpass')
                    downurl = 'https://s1qc.asf.alaska.edu/aux_'+producttype.lower()+'/'+row.filename
                    command = 'wget --user '+ asfuser +' --password '+asfpass+' -O '+outfile+' '+downurl+' 2>/dev/null'
                    rc = os.system(command)
                    #r = requests.get(downurl, allow_redirects=True, auth=HTTPBasicAuth(asfuser, asfpassword))
                    #if r.status_code == 200:
                    #    f = open(outfile, 'wb')
                    #    f.write(r.content)
                    #    f.close()
                    if os.path.exists(outfile):
                        if os.stat(outfile).st_size < 500:
                            os.remove(outfile)
                    if os.path.exists(outfile):
                        print('(probably) ok')
                    else:
                        print('failed also from ASF using wget, sorry')
                        #' - status: '+str(r.status_code))
                except:
                    print('failed also from ASF, sorry')
            os.remove(lockfile)
        if os.path.exists(outfile):
            existing.append(outfile)
    return existing

# get orbit files using eof (they should update once it gets into the new Copernicus cloud... 03/2021)
def updateOrbForZipfile(zipFile, orbdir = os.environ['ORB_DIR']):
    zipFile = os.path.basename(zipFile)
    orbFiles = ''
    #try:
    #    orbFiles = download_eofs(sentinel_file=zipFile, save_dir='.')
    #except:
    #    print('orbits not found using sentineleof. using custom (raw) approach')
    try:
        (startdate, enddate) = strpStrtEndTimes(zipFile)
        startdate = startdate - pd.Timedelta('1 day')
        enddate = enddate + pd.Timedelta('1 day')
        downloadOrbits_CopCloud(startdate, enddate, 'POEORB')
        downloadOrbits_CopCloud(startdate, enddate, 'RESORB')
        return True
    except:
        print('not succeeded')
        return False
    #
    #
    # old unused lines
    if not orbFiles:
        print('orbits not found')
        return False
    for orbF in orbFiles:
        a = os.path.basename(orbF)
        sensor = a.split('_')[0]
        orbType = a.split('_')[3]
        wheresavePath = os.path.join(orbdir, sensor, orbType, a)
        shutil.move(a, wheresavePath)
    #should be anyway just one orbit file...
    #return a
    return wheresavePath


################################################################################
# Get Orbit Url
################################################################################
# 2023/08: seems not working anymore!
def getOrbUrl(sat,prodType,startTime,endTime):
    """Trys to find a valid url for an orbit to download based on the:
    satalite - S1A or S1B
    prodType  (orbit product type) - POEORB or RESORB
    start and end times - note should be date time objects
    """

    #baseUrl = 'https://qc.sentinel1.eo.esa.int/api/v1/?'
    baseUrl = 'https://qc.sentinel1.copernicus.eu/api/v1/?'

    satDict = {'S1A':'S1A',
            'S1B':'S1B'}

    prodDict = {'POEORB':'AUX_POEORB',
            'RESORB':'AUX_RESORB'}

    try:
        prodQry = 'product_type='+prodDict[prodType]
    except KeyError:
        raise invalidSatType(prodType)
    except:
        logger.error("Unexpected error creating product subquery")
        raise

    try:
        satQry = 'sentinel1__mission='+satDict[sat]
    except KeyError:
        raise invalidProductType(sat)
    except:
        logger.error("Unexpected error creating product subquery")
        raise

    #Note the "lt" = less than
    startTimeQry = dt.datetime.strftime(startTime,'validity_start__lt=%Y-%m-%dT%H:%M:%S')
    #Note the "gt" = greater than
    stopTimeQry = dt.datetime.strftime(endTime,'validity_stop__gt=%Y-%m-%dT%H:%M:%S')

    #So we only get the latest orbit
    constrantQry = 'ordering=-creation_date&page_size=1'

    #Combine queries into the query url
    qryUrl = baseUrl+satQry+'&'+prodQry+'&'+startTimeQry+'&'+stopTimeQry+'&'+constrantQry

    
    resp = requests.get(qryUrl)

    resp_json = resp.json()

    if resp_json['count']:
        #Strip out the url if found
        return resp_json['results'][0]['remote_url']
    else:
        raise orbUrlNotFnd(prodType,startTime,endTime)


################################################################################
# Download the orbit file to directory
################################################################################
def downloadOrbit(url,outDir):
    """ Tries to download the orbit in url into the output directory"""    
    eofPat = '.*(S1.*\.EOF)'

    #Strip file name from the url
    mtch = re.search(eofPat,url)

    if mtch:
        outputFile = mtch.groups()[0]
    else:
        raise invalidFileUrl(url)
        
    resp = requests.get(url)

    if resp.ok: #if downloaded write to file
        with open(outDir+'/'+outputFile,'w') as f:
            f.write(resp.content.decode())

        return outputFile
    else:
        raise urlDownloadFailed(url)

################################################################################
# Find orbit file in directory between start and end dates
################################################################################
def findValidOrbFile(baseDir,sat,startTime,endTime):
    """ Searches baseDir for an orbit file which is valid for given satalite,
    start time and end time. Returns None if no file found
    Args:
        baseDir (str): where to search, e.g. /gws/nopw/j04/nceo_geohazards_vol1/orbits.old
        sat (str): either 'S1A' or 'S1B'
        startTime, endTime (dt.datetime)
    """

    timesStrpPat = '.*V(\d{8}T\d*)_(\d{8}T\d*)\.EOF'
    timeStrpPat = '%Y%m%dT%H%M%S'

    satDict = {'S1A':'S1A',
            'S1B':'S1B'}

    try:
        satQry = satDict[sat]
    except KeyError:
        raise invalidProductType(sat)
    except:
        logger.error("Unexpected error creating product subquery")
        raise

    for f in os.listdir(baseDir):
        #Strip out the times in text
        mtch = re.search(satQry+timesStrpPat,f)

        if mtch: # if a match is found (i.e. valid orbit filename)
            #convert to datetime's
            vldTimeStart = dt.datetime.strptime(mtch.groups()[0],timeStrpPat)
            vldTimeStop = dt.datetime.strptime(mtch.groups()[1],timeStrpPat)

            #Check if the orbit time includes the measurement time period
            if (vldTimeStart<=startTime) and (vldTimeStop>=endTime):
                #sometimes the files are wrong and 0 size... so:
                if os.path.getsize(os.path.join(baseDir,f)) < 100:
                    try:
                        os.remove(os.path.join(baseDir,f))
                    except:
                        print('ERROR - the orbit file {} is erroneous but you do not have rights to update it'.format(f))
                else:
                    return f #and return

    return None #otherwise return none

################################################################################
# Get valid orbit file for product
################################################################################
def getValidOrbFile(localOrbDir,prodFile):
    """ Function which searched localOrbDir for a valid orbit file for prodFile.
    If none is found, it will try to download a valid orbit file. It will first 
    look for/download a precise orbit file before falling back on a restituted
    orbit file"""
    
    sat = strpSat(prodFile)
    (startTime,endTime) = strpStrtEndTimes(prodFile)

    #Sub function which searches locally and remotely for a valid orbit

        
    def getOrbFile(orbType):
        baseDir = localOrbDir+'/'+orbType+'/'
        localOrb = findValidOrbFile(baseDir,sat,startTime,endTime)
        if localOrb:
            logger.info("Found local {1} type orbit file - {0}".format(localOrb,orbType))
            return baseDir+localOrb
        else:
            return None
            #logger.info("Couldn't find local {0} type orbit file, locing for url...".format(orbType))
            #orbUrl = getOrbUrl(sat,orbType,startTime,endTime)
            #logger.info("Found URL {0}".format(orbUrl))
            #localOrb = downloadOrbit(orbUrl,baseDir)
            #logger.info("Downloaded {0} type orbit file - {1}".format(orbType,localOrb))
            #return baseDir+localOrb
    
    logger.info("Looking for precise orbit file")
    a = getOrbFile('POEORB')
    if a:
        logger.info('found it')
        return a
    else:
        logger.info("trying to download the files")
        rc =  updateOrbForZipfile(prodFile)
        a = getOrbFile('POEORB')
        if a:
            logger.info('found POEORB')
            return a
        else:
            logger.info("Could not find precise orbit, trying restituted orbit...")
            a = getOrbFile('RESORB')
            if a:
                logger.info('found it')
                return a
            else:
                logger.warning("not downloaded")
                return None

