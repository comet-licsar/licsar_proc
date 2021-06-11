#!/usr/bin/env python
"""
This is a tool intended to be used with LiCS database and can be used to identify
which files are required to build a frame and make nla requests to restore them

Usage:

    LiCSAR_0_getFiles.py -f FRAME -s STARTDATE -e ENDDATE -r -z ZIPFILE.LIST

    -f : use FRAME as frame name
    -s : find files with aquisition date after STARTDATE, in the format YYYY-MM-DD. Default 2016-01-1
    -e : find files with aquisition date before ENDDATE, in the format YYYY-MM-DD. Default curent date
    -r : make nla request, ommit if not wanted
    -z : ZIPFILE.LIST to write output to, ommit if not desired
    -b : batch mode to use with nla request:
            M - by month (default)
            W - by week
    -d : (optional) deletion date for the nla request in the format YYYY-MM-DD
    -n : (optional) notify using user's default email. 
    -N : (optional) notify using given email
"""

################################################################################
#Imports
################################################################################

import LiCSquery as lq
import re
import datetime as dt
import sys
import os.path as path
import nla_client_lib as nla
from getopt import getopt
import pandas as pd
from LiCSAR_lib import s1data

################################################################################
#Exceptions
################################################################################

class badFrameError(Exception):
    def __init__(self,frame):
        self.frame = frame
        self.code = 0
        self.msg = '{0}'
    def __str__(self):
        return str(self.code) + ":" + self.msg.formate(self.frame)

class notOnNLAError(Exception):
    def __init__(self,file):
        self.file
        self.code = -2001
        self.msg = 'Could not find {0} on the NLA tape archive'
    def __str__(self):
        return str(self.code) + ":" + self.msg.formate(self.file)

class undefinedFrameError(badFrameError):
    def __init__(self,frame):
        super(undefinedFrameError,self).__init__(frame)
        self.code = -1001
        self.msg = "{0} is not defined in the LiCS database"

class noFrameGivenError(badFrameError):
    def __init__(self,frame):
        super(noFrameGivenError,self).__init__(frame)
        self.code = -1002
        self.msg = "Frame = {0} - no frame provided"



################################################################################
# Main program to read frame name and request files from nla
################################################################################
def main(argv=None):

    def check_nla(zipFile):
        if not path.exists(zipFile):
            if nla.ls(zipFile,'UDTAR'):
                return True
            else:
                raise notOnNLAError(zipFile)
        return False
    def make_nla_request(fileSeries,userEmail=None):
        if not fileSeries.empty:
            print('getting previous request ids')
            reqInfo = nla.list_requests()
            priorReqs = set([req['id'] for req in reqInfo['requests']])
            print(reqInfo['requests'])
            if not userEmail:
                userEmail = reqInfo['email']
            label = frameName + ': ' + fileSeries.index[0].strftime('%Y-%m-%d') + ' -> ' + fileSeries.index[-1].strftime('%Y-%m-%d')
            print('Requesting: '+label)
            nla.make_request( files = list( fileSeries ),
                    label = label,
                    retention=retentionDate)
            print('nla requested')
            reqInfo = nla.list_requests()
            postReqs = set([req['id'] for req in reqInfo['requests']])
            curReq = postReqs-priorReqs
            if len(curReq) != 0:
                curReq = curReq.pop()
            if curReq:
                print("Created request {0} with label {1}".format( curReq, label))
                if notify:
                    nla.update_request(curReq,notify_first=userEmail,notify_last=userEmail)
                    print("Added notify to requests {0} with email {1}".format(curReq,userEmail))
            else:
                print("Request was not created... Please check manually for label {0}".format(
                        label
                        ))
            return curReq
        else:
            return None

############################## Defaults

    frameName = None
    zipListFile = None
    makeNLAReq = False
    batchMode = 'M'
    startDate = dt.datetime.strptime('2016-01-01','%Y-%m-%d')
    endDate = dt.datetime.now()
    retentionDate = None
    notify = False
    userEmail = None

############################## Parse arguments

    opts,args = getopt(argv[1:],'hrnf:z:b:s:e:N:')

    for o,a in opts:
        if o == '-h':
            print(__doc__)
            return 0
        elif o == '-f':
            frameName = a
        elif o == '-r':
            makeNLAReq = True
        elif o == '-z':
            zipListFile = a
        elif o == '-b':
            batchMode = a
        elif o == '-s':
            startDate = dt.datetime.strptime(a,"%Y-%m-%d")
        elif o == '-e':
            endDate = dt.datetime.strptime(a,"%Y-%m-%d")
        elif o == '-d':
            retentionDate = dt.datetime.strptime(a,"%Y-%m-%d")
        elif o == '-n':
            notify = True
        elif o == '-N':
            notify = True
            userEmail = a
        else:
            print("unexpected argument")
            print(__doc__)
            return 1

############################## Check Frame definition

    if frameName == None:
        raise noFrameGivenError(frameName)
    if not lq.check_frame(frameName):
        raise undefinedFrameError(frameName)

    print("Found frame definition in LiCS database - reading file list")

############################## Check files using scihub -> NLA

    print('checking for S1 data not ingested to licsinfo db')
    s1dataa = s1data.check_and_import_to_licsinfo(frameName, startDate.date(), endDate.date())
    print('check for existing S1 data finished')
    #if not s1dataa:
    #    print('no data to download found, quitting')
    #    return False

############################## Get file list
    print('getting file list')
    frameFilesTable = lq.get_frame_files_period(frameName,startDate.strftime('%Y-%m-%d'),endDate.strftime('%Y-%m-%d'))

    print("Stripping file list down to unique zipfiles")

    acq_dates = [f[1] for f in frameFilesTable]
    files = [f[3] for f in frameFilesTable]
    zipFiles = [re.sub('.metadata_only','',fI) for fI in files]
    
    #fix for zipFiles that are not in /neodc:
    for zipf in zipFiles:
        if 'neodc' not in zipf:
            removeddate = acq_dates.pop(zipFiles.index(zipf))
            zipFiles.remove(zipf)
    
    if s1dataa:
        print('correcting paths for {} missing files'.format(len(s1dataa)))
        for s1f in s1dataa:
            s1neodc = s1data.get_neodc_path_images(s1f.split('.')[0])[0]
            if path.exists(s1neodc.replace('.zip','.manifest')):
                zipFiles.append(s1neodc)
                #files.append(s1neodc)
                s1f_date = dt.datetime.strptime(s1f.split('_')[5].split('T')[0],"%Y%m%d")
                acq_dates.append(s1f_date.date())
    
    filesDF = pd.DataFrame({'files':zipFiles,'onTape':False},index=pd.to_datetime(acq_dates))
    filesDF = filesDF.drop_duplicates()
    
    # this is to correct for duplicate images, see explanation at lq.get_frame_files_period
    pom=''
    pomDF=filesDF
    for index in pomDF.index.unique():
        for file in pomDF.loc[index]['files']:
            #if the previous field had the same base-name, drop this one
            if pom==file[:-9]:
                filesDF=filesDF[filesDF.files != file]
            else: pom=file[:-9]
############################## Write out file list

    if zipListFile:
        print("Writing zip file list to {0}".format(zipListFile))
        with open(zipListFile,'w') as f:
            # for zipFile in zipFiles:
            for date,zipFile in filesDF['files'].items():
                f.write(zipFile+"\n")
    else:
        print("No zip file list filename provided - not writing file list")

############################## Work out which files are on tape and which are not
    print('checking for files existing on Tape')
    filesDF['onTape'] = filesDF['files'].map(check_nla)
    fileSeries = filesDF.loc[filesDF['onTape'],'files']
    fileSeries = fileSeries.sort_index()
    print('There are {0} files to be requested'.format(str(len(fileSeries))))
    existing = len(filesDF['files'])-len(fileSeries)
    if existing > 0:
        print('Seems there are {0} files already existing on disk'.format(str(existing)))
        print('WARNING - we will not update their NLA expiry date')
############################## Make batches of NLA request
    if makeNLAReq:
        print("Making NLA file request")
        requestNumber = make_nla_request(fileSeries,userEmail=userEmail)
        # previously we were separating files to several requests
        # this is not really effective towards NLA system
        # so only one NLA request is now applied
        #print("Making NLA file requests")
        #groupedFiles = filesDF.loc[filesDF['onTape'],'files'].resample(batchMode)
        #requestList = groupedFiles.aggregate(make_nla_request,userEmail=userEmail)

################################################################################
# execute main function
################################################################################
if __name__ == "__main__":
    main(argv=sys.argv)

