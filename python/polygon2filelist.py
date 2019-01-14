#!/usr/bin/env python
"""
Tool which levereges the LiCSAR database to convert a polygon file to lists of 
files stored on tape (NLA) and on disk.

usage:

    polygon2filelist -p polygonfile -d ondisklist -t ontapelist
"""

################################################################################
# imports
################################################################################

import LiCSquery as lq
import re
import os.path as path
import sys
import copy
from getopt import getopt


################################################################################
# Constants
################################################################################
coordRE="(-?\d+\.\d+)\s?,?\s?(-?\d+\.\d+)" # Coord regualar expresion 

################################################################################
# Get polygon from file
################################################################################
def getPolygonFromFile(polyFile):
    
    # Setup polygon set
    polygon = list()

    #iterate over polygon file
    with open(polyFile) as pf:
        for line in pf:
            mtch = re.search(coordRE,line)
            coord = [float(numStr) for numStr in mtch.groups()]
            polygon.append(coord)

    return polygon

################################################################################
# Get polygon long/lat min max
################################################################################
def getPolygonMinMax(polygon):

    polyLon = [coord[0] for coord in polygon]
    polyLat = [coord[1] for coord in polygon]

    return min(polyLon),max(polyLon),min(polyLat),max(polyLat)

################################################################################
# Convert polygon limits to file list
################################################################################
def getFilesFromPolygonMinMax(polyMinMax,track):

    files = set()
    brstTracks = set()

    brstIDs = [brst[0] for brst in lq.get_bursts_in_polygon(*polyMinMax)]

    for brstID in brstIDs:
        curTrack = int(re.search('(\d+)_',brstID).groups()[0])
        if not track:
            brstTracks.add(curTrack)
        brstFiles = [bf[1] for bf in lq.get_files_from_burst(brstID)]
        if curTrack==track:
            files |= set(brstFiles)

    return list(files),list(brstTracks)

################################################################################
# Convert file list to files on disk
################################################################################
def splitFileLists(files):

    filesOnDisk = []
    filesOnTape = []
    zipFiles = set([re.sub('.metadata_only','',fI) 
            for fI in files])
    for zI in zipFiles:
        if path.exists(zI):
            filesOnDisk.append(zI)
        else:
            filesOnTape.append(zI)

    return filesOnDisk,filesOnTape


################################################################################
# Write file list to file
################################################################################
def writeFileList(files,fileListFile):
    with open(fileListFile,'w') as flf:
        for fI in files:
            flf.write(fI+'\n')

################################################################################
# Main program to read polygon file and create a list of zip files
################################################################################
def main(argv=None):

    polygonFile = None
    flsOnDskFile = None
    flsOnTpeFile = None
    burstTrack = None

    opts,args = getopt(argv[1:],'hp:z:b:')

    for o,a in opts:
        if o == '-h':
            print(__doc__)
            return 0
        elif o == '-p':
            polygonFile = a
        elif o == '-z':
            zipListFile = a
        # elif o == '-d':
            # flsOnDskFile = a
        # elif o == '-t':
            # flsOnTpeFile = a
        elif o == '-b':
            burstTrack = int(a)
        else:
            print("unexpected argument")
            print(__doc__)
            return 1

    if len(opts)==0:
        print(__doc__)
        return 0
       
    if (not polygonFile):
        print("polygon file required")
        return 2
    if (not zipListFile):
        print("output file for on disk list required")
        return 3
    
    # if (not flsOnDskFile):
        # print "output file for on disk list required"
        # return 3

    # if (not flsOnTpeFile):
        # print "output file for on tape list required"
        # return 4

    polygon = getPolygonFromFile(polygonFile)
    polygonLim = getPolygonMinMax(polygon)
    if burstTrack:
        files,bts = getFilesFromPolygonMinMax(polygonLim,burstTrack)
        # filesOnDisk,filesOnTape = splitFileLists(files)

        writeFileList(files,zipListFile)
        # writeFileList(filesOnDisk,flsOnDskFile)
        # writeFileList(filesOnTape,flsOnTpeFile)
    else:
        files,bts = getFilesFromPolygonMinMax(polygonLim,burstTrack)
        print("Burst Track required, pick one of the following:")
        for bt in bts:
            print(bt)

################################################################################
# execute main function
################################################################################
if __name__ == "__main__":
    main(argv=sys.argv)

