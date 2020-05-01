import os
import sys
import subprocess as subp
import re

import sys, traceback

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
    open(filee, 'a').close()

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
