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
    with open(fileList) as flf:
        for line in flf:
            dateMatch = re.search('/20\d\d/\d+/\d+/',line)
            dates.append(
                    dt.datetime.strptime(dateMatch.group(),'/%Y/%m/%d/')
                    )
            files.append(line.rstrip())

    return pd.DataFrame({'files':files},index=pd.to_datetime(dates))

################################################################################
# Make batch requests from the NLA
################################################################################
def requestBatchFromNLA(fileList,freq,batchLabel=None,retention=None):
    if not batchLabel:
        batchLabel = 'Batch'
    grouped = fileList.resample(freq)
    grouped.ngroups #dunno - just works
    for indices,group in grouped:
        nla.make_request( files = list( group[ 'files' ] ),
                label = batchLabel + '-' + indices.strftime( '%Y-%m-%d' ),
                retention=retention )

################################################################################
# main function
################################################################################
def main(argv=None):

    fileListFile = None;
    startDate = None;
    endDate = None;
    blockDate = None;
    label = None;
    retention = None;
    frequency = None;

    opts,args = getopt(argv[1:],'hf:s:e:b:l:r:o:')

    for o,a in opts:
        if (o=='-f'):
            fileListFile = a
        elif (o=='-s'):
            startDate = a
        elif (o=='-e'):
            endDate = a
        elif (o=='-b'):
            blockDate = a
        elif (o=='-l'):
            label = a
        elif (o=='-r'):
            retention = dt.strptime(a,'%Y-%m-%d')
        elif (o=='-o'):
            frequency = a
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

    requestBatchFromNLA(fileList,frequency,label,retention)

################################################################################
# execute main function
################################################################################
if __name__ == "__main__":
    main(argv=sys.argv)


