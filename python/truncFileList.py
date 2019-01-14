#!/usr/bin/env python

import nla_client_lib as nla
import datetime as dt
import re
import sys
from getopt import getopt
import pandas as pd

################################################################################
# Load file list
################################################################################
def loadFileList(fileList):
    files = []
    dates = []
    tracks = []
    with open(fileList) as flf:
        for line in flf:
            dateMatch = re.search('/20\d\d/\d+/\d+/',line)
            dates.append(
                    dt.datetime.strptime(dateMatch.group(),'/%Y/%m/%d/')
                    )
            files.append(line.rstrip())

    return pd.DataFrame({'files':files},index=pd.to_datetime(dates))

def writeFileList(fileList,fileListFile):
    files = list(fileList['files'])
    with open(fileListFile,'w') as flf:
        for fI in files:
            flf.write(fI+'\n')

################################################################################
# main function
################################################################################
def main(argv=None):

    fileListFile = None;
    startDate = None;
    endDate = None;
    blockDate = None;
    truncListFile = None;

    opts,args = getopt(argv[1:],'hf:s:e:b:t:')

    for o,a in opts:
        if (o=='-f'):
            fileListFile = a
        elif (o=='-s'):
            startDate = a
        elif (o=='-e'):
            endDate = a
        elif (o=='-b'):
            blockDate = a
        elif (o=='-t'):
            truncListFile = a
        elif (o=='-h'):
            print(__doc__)
            return 0
        else:
            print("unkown option")
            return 1
    
    if not fileListFile:
        print("file list required")
        return 2

    fileList = loadFileList(fileListFile)
    
    if blockDate:
        fileList = fileList[blockDate]
    elif startDate and endDate:
        fileList = fileList[startDate:endDate]
    elif startDate:
        fileList = fileList[startDate:]
    elif endDate:
        fileList = fileList[:endDate]

    writeFileList(fileList,truncListFile)

################################################################################
# execute main function
################################################################################
if __name__ == "__main__":
    main(argv=sys.argv)
