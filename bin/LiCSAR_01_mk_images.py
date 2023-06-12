#!/usr/bin/env python
"""
Identify files related to bursts and frames, and merge into Gamma SLC format.

========
Overview
========
This script is part of the main LiCSAR processing chain. Based on the frame 
and defined date range, it searches the database for related bursts and 
files. For each acquistion date, it unzips the raw files and converts from
SAFE into Gamma format. It then extracts the relevant bursts from each file
and merges them into a new file. Finally it mosaics the subswaths together 
and updates the orbit state vectors. Output logs from individual steps will 
be stored in the log directory in the main processing directory. A high 
level overview of the processing parameters and results is given in the main
processing directory, named '<framename>-mk_images-report.txt'.

=========
Changelog
=========
2020+ various changes (Milan Lazecky) + added cross-pol extraction (for test only)
July 2016: Original implementation (Karsten Spaans, Uni of Leeds)

=====
Usage
=====
LiCSAR_01_mk_images.py -n -i <ignorelist> -f <framename> -d </path/to/processing/location> -s <startdate> -e <enddate> -m <masterdate> -a <azlooks> -r <rglooks>

    -f    Name of the frame to be processed in database
    -d    Path to the processing location. 
    -s    Earliest date (YYYYMMDD) of images that should be processed. 
          If not given, all images before end date are processed.
    -e    Latest date (YYYYMMDD) of images that should be processed. 
          If not given, all images after start date are processed. 
          If neither start nor end date are given, all images are processed.
    -a    Integer azimuth multilook factor, defaults to 4
    -r    Integer range multilook factor, defaults to 20
    -m    Date of master image. If not given, it will give a list of 
          available dates
    -p    Polygon file containing longitude and latitude coordinates which
          define the frame.
    -z    Zip file list containing a list of zip files of the raw data.
    -y    Batch mode flag, set to 1 if running from terminal. In Batch mode 
          polygon and Zip file list must be provided!
    -t    Track number
    -n    no master - extract without master date
    -b    File containing burst ids
    -j    Autoproccessing job id
    -i    load ignore dates file. Will ignore dates in this file (YYYYMMDD) can be generated using ls SLC >> ignore.list
    -l <file> only extract specific dates in file. Date format is (YYYYMMDD)
    -T <file> Report (text) to write to. Defaults to FRAME-mk_images-report.txt
    -c    Clean mosaics - This will remove mosaicked slc's, leaving only the subswath slcs.
            This will reduce disk usage during processing.
    -x    Run in test crosspol mode (i.e. extract only VH or HV data)

"""

################################################################################
# Import relevent packages
################################################################################
import getopt
import os
import sys
import numpy as np
import datetime as dt
# dir_path = os.path.dirname(os.path.realpath(__file__)) # This should be gone...
# sys.path.append(dir_path[:-4]+'/lib')
# sys.path.append(dir_path[:-4]+'/LiCSdb')
import global_config as gc
from LiCSAR_lib.mk_imag_lib import *
from LiCSAR_lib.LiCSAR_misc import *

################################################################################
#Main program function
################################################################################
def main(argv=None):
    if argv == None:
        argv = sys.argv
    print(argv)

############################################################ Param init.
    framename=[]
    burstidfile = []
    polygonfile = []
    ziplistfile = []
    procdir=[]
    startdate=[]
    enddate=[]
    masterdate=[]
    ignoreListFile=None
    dateListFile=None
    trackno=[]
    job_id = -1
    skipMstr = False
    reportfile = None
    removeSlcs = False
    make_local_config = False
    test_crosspol = False

############################################################ Parse program args.
    try:
        try:
            opts, args = getopt.getopt(argv[1:], 
                    "vhnxci:f:d:s:e:m:p:b:o:t:j:a:r:z:y:l:T:",
    		    ["version", "help"])
        except getopt.error as msg:
            raise Usage(msg)
        for o, a in opts:
            if o == '-h' or o == '--help':
                print(__doc__)
                return 0
            elif o == '-v' or o == '--version':
                print("")
                print("Current version: %s" % gc.config['VERSION'])
                print("")
                return 0
            elif o == '-f':
                framename = a
            elif o == '-d':
                procdir = a
            elif o == '-s':
                startdate = dt.date(int(a[:4]),int(a[4:6]),int(a[6:8]))
            elif o == '-e':
                enddate = dt.date(int(a[:4]),int(a[4:6]),int(a[6:8]))
            elif o == '-m':
                masterdate = dt.date(int(a[:4]),int(a[4:6]),int(a[6:8]))
            elif o == '-p':
                polygonfile = a
            elif o == '-b':
                burstidfile = a
            elif o == '-t':
                tracknumber = a
            elif o == '-n':
                skipMstr = True
            elif o == '-j':
                job_id = int(a)
            elif o == '-a':
                if int(gc.azlks) != int(a):
                    make_local_config = True
                gc.azlks = int(a)
            elif o == '-r':
                if int(gc.rglks) != int(a):
                    make_local_config = True
                gc.rglks = int(a)
            elif o == '-z':
                ziplistfile = a
            elif o == '-y':
                gc.batchflag = bool(a)
            elif o == '-i':
                ignoreListFile = a
            elif o == '-l':
                dateListFile = a
            elif o == '-T':
                reportfile = a
            elif o == '-c':
                removeSlcs = True
            elif o == '-x':
                test_crosspol = True
                print('WARNING, extraction of cross-pol data not tested well')
############################################################ Check 4 Crit. Params.
        if not (framename or polygonfile or burstidfile):
            raise Usage('No frame, polygon or set of burst ids given,'\
                    'please use either the -f, -p or -b option!')

        if not procdir:
            raise Usage('No output data directory given, -d is not optional!')

############################################################ Job & batch processing...
        if job_id == -1:
            print("This processing is not outputting any products to the "\
                    "database.")
        
        global lq 
        if gc.batchflag:
            print("Batch processing mode in use --> no database interaction.")
            if (not polygonfile) :
                raise Usage("Polygon file required for batch mode. Use -p "\
                        "option.")
            if (not ziplistfile):
                raise Usage("Ziplist file required for batch mode. Use -z "\
                        "option.")
            import LiCSquery_batch
            lq = LiCSquery_batch.dbquery(ziplistfile, polygonfile)
        else:
            import LiCSquery as lq

    
############################################################ Incorrect usage 
                                                           # err. handling
    except Usage as err:
        print("\nERROR:", file=sys.stderr)
        print("  "+str(err.msg), file=sys.stderr)
        print("\nFor help, use -h or --help.\n", file=sys.stderr)
        return 2
    
############################################################ D.B. Connection Check
    #Check if a DB connection can be established
    if not lq.connection_established():
        print("\nERROR:", file=sys.stderr)
        print("Could not establish a stable database "\
                "connection. No processing can happen.", file=sys.stderr)

        return 1

############################################################ Create Proc. Dir.
    # Check if proc dir can be created
    if not os.path.exists(procdir):
        os.mkdir(procdir)
    if make_local_config:
        lc_file = os.path.join(procdir, 'local_config.py')
        with open(lc_file,'w') as lc:
            lc.write('azlks = {0}\n'.format(str(gc.azlks)))
            lc.write('rglks = {0}\n'.format(str(gc.rglks)))
            #assuming higher resolution, lets create hires geo files
            lc.write('outres = 0.0001\n')

############################################################ Date init.
    # Start and end date needed for database query
    if not startdate:
        startdate=dt.date(2014,10,1)
    if not enddate:
        enddate=dt.date.today()

############################################################ Check crit. params. exist.
    # Check if frame in database
    if framename:
        testlist = lq.check_frame(framename)
        if not testlist:
            print("\nERROR:", file=sys.stderr)
            print("Could not find frame {0} in "\
                    "database".format(framename), file=sys.stderr)
            return 1
    if polygonfile:
        if not os.path.exists(polygonfile):
            print("\nERROR:", file=sys.stderr)
            print("Could not find the given polygon file"\
                    "{0}".format(polygonfile), file=sys.stderr)
            return 1
        if (not trackno) and (not framename) and (not gc.batchflag):
            print("\nERROR:", file=sys.stderr)
            print("No track number given with polygonfile, "\
                    "use -t option to define".format(polygonfile), file=sys.stderr)
            return 1
    if burstidfile:
        if not os.path.exists(burstidfile):
            print("\nERROR:", file=sys.stderr)
            print("Could not find the given burst id file "\
                    "{0}".format(burstidfile), file=sys.stderr)
            return 1
        if (not trackno) and (not framename) and (not gc.batchflag):
            print("\nERROR:", file=sys.stderr)
            print("No track number given with burst id file, "\
                    "use -t option to define".format(polygonfile), file=sys.stderr)
            return 1
    if gc.batchflag:
            pass
    print('\nInput parameters seem okay! Checking frame, bursts and "\
            "files information...')
    



############################################################ load master data from poly./frame
    if framename:
        res = check_bursts(framename,startdate,enddate,lq)
        if res != 1:
            burstlist, filelist, dates = res
        else:
            return 1
    elif polygonfile:
        with open(polygonfile) as pf:
            polygon = np.loadtxt(pf,dtype=np.float32)
            if len(polygon) <4 or len(polygon) > 5 and not(gc.batchflag):
                print("\nERROR:", file=sys.stderr)
                print("Given polygon is not rectangular. "\
                        "Please provide a rectangular polygon, either closed or open.", file=sys.stderr)
                return 1
    elif burstidfile:
        pass
    
    if gc.batchflag:
        # Needs to define:
        # - framename (for labelling etc)
        # - burstlist
        # - zipfilelist
        framename = ''.join(os.path.split(polygonfile)[-1].split('.')[0:-1])
        print('\nDefined framename: {0}'.format(framename))
        burstlist = lq.get_bursts_in_polygon()
        dates = lq.get_dates()
        
    # Load usr specified dates - then find common dates in zipfile list
    if dateListFile:
        usrDates = list()
        with open(dateListFile) as f:
            for line in f:
                usrDates.append(dt.datetime.strptime(line.strip(),'%Y%m%d').date())
        usrDates = set(usrDates)
        if usrDates-set(dates):
            print("Warning following date were not in zipfile list {0}".format(usrDates-set(dates)))
        dates = list(usrDates&set(dates))

    #Dates to ignore....
    if ignoreListFile:
        with open(ignoreListFile) as f:
            for line in f:
                ignrDate = dt.datetime.strptime(line.strip(),'%Y%m%d').date()
                dates.remove(ignrDate)

    
############################################################ Check Master Bursts
    # Check if master is in date list, and has all bursts. Missing bursts ok
    # for slaves, but not for master.
    if not skipMstr:
        if masterdate:
            rc = check_master_bursts(framename,burstlist,masterdate,dates,lq)
            if rc != 0:
                return 1
        else:
            print('\nNo master date given. Please use the -m "\
                    "option to define one of these choices for the master:\n{0}'.format(
                            ', '.join([m.strftime('%Y%m%d') for m in sorted(list(dates))
                                if m != masterdate])), file=sys.stderr)
            return 1
    else:
        print("\nSkipping master burst check")


############################################################ Report File init.
    # Initialize report file
    if not reportfile:
        reportfile = os.path.join(procdir,'{0}-mk_images-report.txt'.format(framename))
    starttime = dt.datetime.now()
    with open(reportfile,'w') as f:
        f.write('LiCSAR_01_mk_images processing started {0}.\n'.format(
            dt.datetime.now().strftime('%Y-%m-%d %H:%M')))
        if skipMstr:
            f.write('\nNo master date')
        else:
            f.write('\nMaster date: {0}\n'.format(masterdate))
        f.write('Start date: {0}\n'.format(startdate))
        f.write('End date: {0}\n'.format(enddate))
        f.write('Range looks: {0}\n'.format(gc.rglks))
        f.write('Azimuth looks: {0}\n'.format(gc.azlks))
        f.write('\n{0} images were found.\n'.format(len(dates)))

############################################################ Date Loop
    # Looping through files
    okcount = 0
    print('\nUnzipping, converting and merging image files per date...')
    for date in dates:
        print('\nProcessing acquisition date {0}:'.format(date))
        #print(polygon)
        #print(polygon)
        if gc.batchflag:
            imburstlist = lq.get_bursts_in_polygon()
            #print(imburstlist)
        elif framename:
            imburstlist = lq.get_frame_bursts_on_date(framename,date)
        #elif polygon:
        #    imburstlist = lq.get_bursts_in_polygon(polygon[:,0].min(),
        #                                           polygon[:,0].max(),
        #                                           polygon[:,1].min(),
        #                                           polygon[:,1].max())
        #    print('image burst list:')
        #    print(imburstlist)
        #else:
        #    imburstlist = []
        #    print('for debug purposes - skipping check for ')
############################################################ if no missing bursts -> Make Frames
        if imburstlist:
            # Check if all bursts of frame are in current image
            missingbursts = [b for b in burstlist if not b in imburstlist]
            if len(missingbursts) < 1:
                # All bursts there, no problem
                print("All bursts for frame {0} seem to be have been acquired "\
                        "on {1}...".format(framename,date))
                rc = make_frame_image(date,framename,imburstlist,procdir, lq, job_id, test_crosspol = test_crosspol) # Key processing function
                if removeSlcs and rc == 0 :
                    slc = os.path.join(procdir,'SLC',date.strftime('%Y%m%d'),
                                        date.strftime('%Y%m%d.slc'))
                    print("Removing mosaiced image {0}".format(slc))
                    os.remove(slc)
                        

############################################################ Deal with missing bursts...
            else:
                # Missing one or more bursts, checking if in problematic loc
                print("One of more  bursts for frame {0} have not been acquired "\
                        "on the date {1}. Missing bursts: {2}".format(framename,
                                date,''.join(['\n'+m[0] for m in missingbursts])))
                print("Checking where missing bursts are.")
                if check_missing_bursts(burstlist,missingbursts):
                    # Missing bursts are at edges, not a problem
############################################################ OK -> Make Frames
                    print("Missing bursts are not in critical location, continuing "\
                            "processing image {0}.".format(date))
                    rc = make_frame_image(date,framename,imburstlist,procdir, lq,job_id, test_crosspol = test_crosspol) # Key processing function
                    if removeSlcs and rc == 0:
                        slc = os.path.join(procdir,'SLC',date.strftime('%Y%m%d'),
                                            date.strftime('%Y%m%d.slc'))
                        print("Removing mosaiced image {0}".format(slc))
                        os.remove(slc)
                else:
                    # Missing bursts in the middle, we're stuffed for 
                    # this date, skip to next
                    print("Missing bursts are in critical location, continuing "\
                            "processing with next image...")
                    rc = 4
        # datedir = os.path.join(procdir,'SLC',date.strftime('%Y%m%d')) # is this needed?
        else:
            print('WARNING: list of bursts is empty - probably not using LiCSInfo?')
            print('trying to go on, no warranties')
            rc = make_frame_image(date,framename,imburstlist,procdir, lq,job_id, test_crosspol = test_crosspol) # Key processing function
        
############################################################ Processing status > report file
        with open(reportfile,'a') as f:
            if rc == 0:
                okcount +=1
            elif rc == 1:
                f.write('\nAcquisition {0} had a problem during the file reading.'.format(date))
#                shutil.rmtree(datedir)
            elif rc == 2:
                f.write('\nAcquisition {0} had a problem during the copying and merging of bursts.'.format(date))
#                shutil.rmtree(datedir)                
            elif rc == 3:
                f.write('\nAcquisition {0} had a problem during the mosaicing of subswaths.'.format(date))
#                shutil.rmtree(datedir)
            elif rc == 4:
                f.write('\nAcquisition {0} has missing bursts in critical location.'.format(date))
#                if os.path.exists(datedir):
#                    shutil.rmtree(datedir)
    try:
        lq.close_db_and_tunnel()
    except:
        print('')
    
############################################################ Finish Report File
    # Finalize report file
    with open(reportfile,'a') as f:
        f.write('\n\nDone. Successfully processed {0} out of {1} images\n'.format(okcount,len(dates)))
    
if __name__ == "__main__":
    sys.exit(main())

