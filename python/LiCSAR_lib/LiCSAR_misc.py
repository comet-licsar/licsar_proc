import os
import sys
import subprocess as subp
import re
import zipfile
import sys, traceback
from osgeo import gdal, gdalconst

class nostdout(object):
    def __enter__(self):
        self.stdout = sys.stdout
        sys.stdout = self
    def __exit__(self, type, value, traceback):
        sys.stdout = self.stdout
        #if type is not None:
        #    # Do normal exception handling
    def write(self, x): pass
#usage:
#with nostdout():
#    DoMyFunction(*args,**kwargs)


################################################################################
#CD class
################################################################################
class cd:
    """Context manager for changing the current working directory, used with
    while to provide a temporary path change"""
    def __init__(self, newPath):
        self.newPath = os.path.expanduser(newPath)

    def __enter__(self):
        self.savedPath = os.getcwd()
        os.chdir(self.newPath)

    def __exit__(self, etype, value, traceback):
        os.chdir(self.savedPath)

################################################################################
# Usage exception class
################################################################################
class Usage(Exception):
    """Usage context manager"""
    def __init__(self, msg):
        self.msg = msg

################################################################################
# Grep function (grep wrapper?)
################################################################################
def grep(arg,file):
    res = subp.check_output(['grep',arg,file])
    return res


def touch(filee):
    if os.path.exists(os.path.dirname(filee)):
        open(filee, 'a').close()
    else:
        print('ERROR - can touch files only in existing folder')


################################################################################
# Grep for only first occurrence (better for text files in python3+)
################################################################################
def grep1(arg,filename):
    file = open(filename, "r")
    res=''
    for line in file:
        if re.search(arg, line):
            res=line
            break
    file.close()
    return res

def grep1line(arg,filename):
    file = open(filename, "r")
    res=''
    for line in file:
        if re.search(arg, line):
            res=line
            break
    file.close()
    if res:
        res = res.split('\n')[0]
    return res

def getipf(parfile):
    line = grep1('IPF', parfile)
    res = float(line.split('IPF')[1].split(')')[0])
    return res

#########################
## sed helper

def sed_rmlinematch(oldstr, infile, dryrun=False):
    '''
    Sed-like line deletion function based on given string..
    Usage: pysed.rmlinematch(<Unwanted string>, <Text File>)
    Example: pysed.rmlinematch('xyz', '/path/to/file.txt')
    Example 'DRYRUN': pysed.rmlinematch('xyz', '/path/to/file.txt', dryrun=True) #This will dump the output to STDOUT instead of changing the input file.
    '''
    linelist = []
    with open(infile) as f:
        for item in f:
            rmitem = re.match(r'.*{}'.format(oldstr), item)
            if type(rmitem) == type(None): linelist.append(item)
    if dryrun == False:
        with open(infile, "w") as f:
            f.truncate()
            for line in linelist: f.writelines(line)
    elif dryrun == True:
        for line in linelist: print(line, end='')
    else:
        print("Unknown option specified to 'dryrun' argument, Usage: dryrun=<True|False>.")
        return False


def sed_replace(oldstr, newstr, infile, dryrun=False):
    '''
    Sed-like Replace function..
    Usage: pysed.replace(<Old string>, <Replacement String>, <Text File>)
    Example: pysed.replace('xyz', 'XYZ', '/path/to/file.txt')
    Example 'DRYRUN': pysed.replace('xyz', 'XYZ', '/path/to/file.txt', dryrun=True) #This will dump the output to STDOUT instead of changing the input file.
    '''
    linelist = []
    with open(infile) as f:
        for item in f:
            newitem = re.sub(oldstr, newstr, item)
            linelist.append(newitem)
    if dryrun == False:
        with open(infile, "w") as f:
            f.truncate()
            for line in linelist: f.writelines(line)
    elif dryrun == True:
        for line in linelist: print(line, end='')
    else:
        print("Unknown option specified to 'dryrun' argument, Usage: dryrun=<True|False>.")
        return False


def is_good_zip(zipf):
    try:
        the_zip_file = zipfile.ZipFile(zipf)
    except:
        print('error opening zip file - marking as bad')
        return False
    try:
        ret = the_zip_file.testzip()
    except:
        print('error testing zip file - marking as bad')
        return False
    if ret is not None:
        print('the file is corrupted')
        #print "First bad file in zip: %s" % ret
        return False
    else:
        return True


################################################################################
# Get information if the file exists or is non-zero
################################################################################
def is_non_zero_file(fpath): 
    return os.path.isfile(fpath) and os.path.getsize(fpath) > 0


# some geometry functions... may or may not be part of LiCSquery?
def get_centre_from_latlon(latlon):
    #latlon should be e.g. [(15.088, -24.56),(14.099, -25.66689)] 
    lats = []
    lons = []
    for coord in range(len(latlon)):
        lats.append(latlon[coord][0])
        lons.append(latlon[coord][1])
    centre = ((min(lats)+max(lats))/2, (min(lons)+max(lons))/2)
    return centre

def get_colat10(lat):
    lat = round(lat * 100)
    if lat < 0:
        lat = abs(lat)
        colat = lat + 90*100
    else:
        colat = 90*100 - lat
    return colat


def reproject_to_match(src_filename, match_filename, dst_filename):
    src = gdal.Open(src_filename, gdalconst.GA_ReadOnly)
    src_proj = src.GetProjection()
    src_geotrans = src.GetGeoTransform()
    srcdatatype = src.GetRasterBand(1).DataType
    #srcdatatype = gdal.GetDataTypeName(srcdatatype)
    nodatav = src.GetRasterBand(1).GetNoDataValue()
    
    # We want a section of source that matches this:
    match_ds = gdal.Open(match_filename, gdalconst.GA_ReadOnly)
    match_proj = match_ds.GetProjection()
    match_geotrans = match_ds.GetGeoTransform()
    wide = match_ds.RasterXSize
    high = match_ds.RasterYSize
    
    # Output / destination
    dst = gdal.GetDriverByName('GTiff').Create(dst_filename, wide, high, 1, srcdatatype) # not working: ['COMPRESS=LZW'])
    dst.SetGeoTransform( match_geotrans )
    dst.SetProjection( match_proj)
    band = dst.GetRasterBand(1)
    band.SetNoDataValue(nodatav)
    
    # Do the work
    gdal.ReprojectImage(src, dst, src_proj, match_proj, gdalconst.GRA_NearestNeighbour)
    del dst # Flush
    return dst_filename

