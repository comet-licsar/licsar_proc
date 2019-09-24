
################################################################################
# Imports
################################################################################

import requests
import re
import datetime as dt
import os
import logging


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

################################################################################
# Get Orbit Url
################################################################################
def getOrbUrl(sat,prodType,startTime,endTime):
    """Trys to find a valid url for an orbit to download based on the:
    satalite - S1A or S1B
    prodType  (orbit product type) - POEORB or RESORB
    start and end times - note should be date time objects
    """

    baseUrl = 'https://qc.sentinel1.eo.esa.int/api/v1/?'

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
            f.write(resp.content)

        return outputFile
    else:
        raise urlDownloadFailed(url)

################################################################################
# Find orbit file in directory between start and end dates
################################################################################
def findValidOrbFile(baseDir,sat,startTime,endTime):
    """ Searches baseDir for an orbit file which is valid for given satalite,
    start time and end time. Returns None if no file found """

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
            logger.info("Couldn't find local {0} type orbit file, locing for url...".format(orbType))
            orbUrl = getOrbUrl(sat,orbType,startTime,endTime)
            logger.info("Found URL {0}".format(orbUrl))
            localOrb = downloadOrbit(orbUrl,baseDir)
            logger.info("Downloaded {0} type orbit file - {1}".format(orbType,localOrb))
            return baseDir+localOrb

    try:
        logger.info("Looking for precise orbit file")
        return getOrbFile('POEORB')
    except orbUrlNotFnd:
        logger.info("Could not find precise orbit, trying restituted orbit...")
        try:
            return getOrbFile('RESORB')
        except orbUrlNotFnd:
            logger.warning("could not find restituted orbit either.")
            raise orbUrlNotFnd
