#!/usr/bin/env python
"""
Unwrap interferograms

========
Overview
========
This script is part of the main LiCSAR processing chain. It requires LiCSAR_03_mk_ifgs.py to have been completed successfully for all slaves. This script filters the wrapped interferograms and performs the unwrapping

=========
Changelog
=========
May 2017: Original implementation (Karsten Spaans, Uni of Leeds)
2020: added parallelisation etc. (M. Lazecky, Uni of Leeds)

=====
Usage
=====
LiCSAR_04_unwrap.py -d </path/to/processing/location>

    -d    Path to the processing location. 
    -j    job id
    -p    polygon file
    -z    zipfile list
    -y    batch mode
    -n    if set to number of processors (>1), it will attempt to parallelise
    -l <file> file containing a list of interferograms to unwrap in format YYYYMMDD_YYYYMMDD
    -T <file> report file to write to. Defaults to report-unwrap.txt
"""

################################################################################
#Imports
################################################################################
import getopt
import os
import sys
import re
import datetime as dt
from gamma_functions import *
import global_config as gc
from LiCSAR_lib.LiCSAR_misc import *
from LiCSAR_lib.unwrp_lib import *

################################################################################
#main function
################################################################################
def main(argv=None):
    if argv == None:
        argv = sys.argv
    
    parallelise = False
    procdir = []
    framename = []
    ifgListFile = None
    job_id = -1
    reportfile = None
    try:
        try:
            opts, args = getopt.getopt(argv[1:], "vhl:f:d:j:n:y:z:p:T:", ["version", "help"])
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
            elif o == '-n':
                num_processors = int(a)
                parallelise = True
            elif o == '-j':
                job_id = int(a)
            elif o == '-y':
                gc.batchflag = bool(a)
            elif o == '-z':
                ziplistfile = a
            elif o == '-p':
                polygonfile = a
            elif o == '-l':
                ifgListFile = a
            elif o == '-T':
                reportfile = a

        if parallelise:
            try:
                from pathos.multiprocessing import ProcessingPool as Pool
                #from multiprocessing import Pool
            except:
                print('error establishing pathos multiprocessing pool - will work without parallelisation')
                parallelise=False
        #if not framename:
        #    raise Usage('No frame name given, -f option is not optional!')
        if not procdir:
            raise Usage('No output data directory given, -d is not optional!')
        if job_id == -1:
            print("This processing is not outputting any products to the database.")
        if gc.batchflag:
            if not(polygonfile) or not(ziplistfile):
                raise Usage("Polygon file AND ziplist file are required in batch mode.")
            import LiCSquery_batch
            global lq
            lq = LiCSquery_batch.dbquery(ziplistfile, polygonfile)
        else:
            import LiCSquery as lq

    except Usage as err:
        print("\nERROR:", file=sys.stderr)
        print("  "+str(err.msg), file=sys.stderr)
        print("\nFor help, use -h or --help.\n", file=sys.stderr)
        return 2

    #Check if a DB connection can be established
    if not lq.connection_established():
        print("\nERROR:", file=sys.stderr)
        print("Could not establish a stable database connection. No processing can happen.", file=sys.stderr)

        return 1

    #if not lq.check_frame(framename):
    #    print("\nERROR:", file=sys.stderr)
    #    print("Frame {0} was not found in database.".format(framename), file=sys.stderr)
    #    return 1

    if not os.path.exists(procdir):
        print("\nERROR:", file=sys.stderr)
        print("Processing directory {0} does not seem to exist.".format(procdir), file=sys.stderr)
        return 1

    ifgdir = os.path.join(procdir,'IFG')
    if not os.path.exists(ifgdir):
        print("\nERROR:", file=sys.stderr)
        print("Can't find IFG directory in processing directory {0}.".format(procdir), file=sys.stderr)
        return 1

    if not os.path.exists(os.path.join(procdir,'log')):
        os.mkdir(os.path.join(procdir,'log'))

    masterdatestr = []
    geodir = os.path.join(procdir,'geo')
    geodirlist = os.listdir(geodir)
    for l in geodirlist:
        if l[-4:] == '.hgt' and l[0] == '2':
            masterdatestr = l[:8]
            
    if masterdatestr:
        masterdate = dt.date(int(masterdatestr[:4]),int(masterdatestr[4:6]),int(masterdatestr[6:8]))
    else:
        print("\nERROR:", file=sys.stderr)
        print("Could not find master date in geocoding directory {0}.".format(geodir), file=sys.stderr)
        return 1

    if not reportfile:
        reportfile = os.path.join(procdir,'report-unwrap.txt')
    starttime = dt.datetime.now()
    with open(reportfile,'w') as f:
        f.write('LiCSAR_04_unwrap processing started {0}.\n'.format(dt.datetime.now().strftime('%Y-%m-%d %H:%M')))
        f.write('\n'.format(masterdate))

    if ifgListFile:
        ifglist = list()
        with open(ifgListFile,'r') as f:
            for line in f:
                mtch = re.search('\d{8}_\d{8}',line.strip())
                if mtch:
                    ifglist.append(mtch.group())
    else:
        ifglist = os.listdir(ifgdir)
    success = 0
    if job_id == -1:
        #print('lq not needed now, closing')
        try:
            lq.close_db_and_tunnel()
        except:
            print('')
    def unwrap_this_ifg(ifg, write_output=False):
        success = 0
        rc = do_unwrapping(masterdatestr,ifg,ifgdir,procdir,lq,job_id)
        if write_output:
            with open(reportfile,'a') as f:
                if rc == 0:
                    success +=1
                if rc ==1:
                    # Filtering
                    f.write('\Interferogram {0} had a problem during the filtering.'.format(ifg)) 
                elif rc == 2:
                    # Unwrapping
                    f.write('\Interferogram {0} had a problem during the unwrapping.'.format(ifg)) 
        #return success
    if parallelise:
        try:
            p = Pool(num_processors)
            p.map(unwrap_this_ifg, ifglist)
            #p.close()
            #p.join()
        except:
            print('some error in parallelisation - trying without it')
            for ifg in ifglist:
                unwrap_this_ifg(ifg, True)
    else:
        for ifg in ifglist:
            unwrap_this_ifg(ifg, True)
    with open(reportfile,'a') as f:
        print('Elapsed time {0}'.format(dt.datetime.now() - starttime))
        f.write('Unwrapping job completed')
    if not job_id == -1:
        try:
            lq.close_db_and_tunnel()
        except:
            print('')


if __name__ == "__main__":
    sys.exit(main())

