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
# Load file list
################################################################################
def loadFileList(fileList):
    files = []
    with open(fileList) as flf:
        for line in flf:
            files.append(line.rstrip())

    return files

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

    zipListFile = None
    flsOnDskFile = None
    flsOnTpeFile = None
    burstTrack = None

    opts,args = getopt(argv[1:],'hz:t:d:')

    for o,a in opts:
        if o == '-h':
            print(__doc__)
            return 0
        elif o == '-z':
            zipListFile = a
        elif o == '-d':
            flsOnDskFile = a
        elif o == '-t':
            flsOnTpeFile = a
        else:
            print("unexpected argument")
            print(__doc__)
            return 1

    if len(opts)==0:
        print(__doc__)
        return 0
       
    if (not zipListFile):
        print("output file for on disk list required")
        return 2
    
    if (not flsOnDskFile):
        print("output file for on disk list required")
        return 3

    if (not flsOnTpeFile):
        print("output file for on tape list required")
        return 4
    files = loadFileList(zipListFile)

    filesOnDisk,filesOnTape = splitFileLists(files)

    writeFileList(filesOnDisk,flsOnDskFile)
    writeFileList(filesOnTape,flsOnTpeFile)

################################################################################
# execute main function
################################################################################
if __name__ == "__main__":
    main(argv=sys.argv)

