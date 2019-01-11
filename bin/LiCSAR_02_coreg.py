#!/usr/bin/env python
"""
Coregister slave images to master

========
Overview
========
This script is part of the main LiCSAR processing chain. It requires step LiCSAR_01_mk_images to have been run successfully on all slaves and the master, and step LiCSAR_01_mk_crop_extDEM to have been finished. This script geocodes the master using the DEM, and coregisters and resamples each slave to the master. A high level overview of the processing parameters and results is given in the main processing directory, named '<framename>-coreg-report.txt'.

=========
Changelog
=========
August 2016: Original implementation (Karsten Spaans, Uni of Leeds)


=====
Usage
=====
LiCSAR_02_coreg.py -f <framename> -d </path/to/processing/location> -m <masterdate>

    -f    Name of the frame to be processed in database
    -d    Path to the processing location. 
    -m    Date of master image. If not given, it will give a list of available dates
    -b    =1 for Batch mode,  =0 (default) for automated (depriciated)
    -y    batch mode flag
    -i    skip geocoding of master
    -j    job id
    -z    zip file list
    -p    polygon file to use
    -l <file>   File containing a list of dates to coregester in the formaty (YYYYMMDD)
    -k    Keep SLCs
    -C    Delete mosaiced rslc's. This will also reduce the runtime space usuage, but they
          will have to be recreated from subswathes (so don't use this option with -c)
    -c    Clean subswathe rslc's. This will reduce the final disk usuage, but will prevent 
          these dates being used for spectral diversity.
    -T <file> report file to use. Defaults to FRAME-coreg-report.txt
    -R    Force clean recreation of rslcs - by default if it finds an existing lookup table
          it will use this instead

"""

################################################################################
#Import Modules
################################################################################
import getopt
import os
import shutil
import sys
import numpy as np
import datetime as dt
from glob import glob
# import subprocess as subp

#Import Gamma functions
# from gamma_functions import *
#And import the global configuration
import global_config as gc
#Python debugger module
# import pdb
from LiCSAR_lib.LiCSAR_misc import *
from LiCSAR_lib.coreg_lib import *

################################################################################
# Main program function
################################################################################
def main(argv=None):
############################################################ Default parameters
    if argv == None:
        argv = sys.argv
    procdir = []
    framename = []
    ignrMstr = False
    batchflag = 0
    job_id = -1
    coregListFile = None
    keepSLCs = False
    removeMosaic = False
    removeSSRslc = False
    reportfile = None
    forceRecreate = False
############################################################ Parse arguments
    try:
        try:
            opts, args = getopt.getopt(argv[1:], "vhkiCcRl:f:d:m:j:y:z:p:T:", ["version", "help"])
        except getopt.error, msg:
            raise Usage(msg)
        for o, a in opts:
            if o == '-h' or o == '--help':
                print __doc__
                return 0
            elif o == '-v' or o == '--version':
                print ""
                print "Current version: %s" % gc.config['VERSION']
                print ""
                return 0
            elif o == '-f':
                framename = a
            elif o == '-d':
                procdir = a
            elif o == '-i':
                ignrMstr = True
            elif o == '-m':
                masterdate = dt.date(int(a[:4]),int(a[4:6]),int(a[6:8]))
            elif o == '-j':
                job_id = int(a)
            elif o == '-y':
                gc.batchflag = int(a)
            elif o == '-z':
                ziplistfile = a
            elif o == '-p':
                polygonfile = a
            elif o == '-l':
                coregListFile = a
            elif o == '-k':
                keepSLCs = True
            elif o == '-C':
                removeMosaic = True
            elif o == '-c':
                removeSSRslc = True
            elif o == '-T':
                reportfile = a
            elif o == '-R':
                forceRecreate = True
        
############################################################ Ensure critical parameters are present
        if not framename:
            raise Usage('No frame name given, -f option is not optional!')
        if not procdir:
            raise Usage('No output data directory given, -d is not optional!')
        if job_id == -1:
            print "This processing is not outputting any products to the database."
############################################################ Setup processing database
        if gc.batchflag:
            if not(polygonfile) or not(ziplistfile):
                raise Usage("Polygon file AND ziplist file are required in batch mode.")
            import LiCSquery_batch
            global lq
            lq = LiCSquery_batch.dbquery(ziplistfile, polygonfile)
        else:
            import LiCSquery as lq
            

    except Usage, err:
        print >>sys.stderr, "\nERROR:"
        print >>sys.stderr, "  "+str(err.msg)
        print >>sys.stderr, "\nFor help, use -h or --help.\n"
        return 2

############################################################ Ensure we have a connection to the database
    #Check if a DB connection can be established
    if not lq.connection_established():
        print >> sys.stderr, "\nERROR:"
        print >> sys.stderr, "Could not establish a stable database connection. No processing can happen."

        return 1

############################################################ Check that parameters are in database
    if not lq.check_frame(framename):
        print >>sys.stderr, "\nERROR:"
        print >>sys.stderr, "Frame {0} was not found in database.".format(framename)

    if not os.path.exists(procdir):
        print >>sys.stderr, "\nERROR:"
        print >>sys.stderr, "Processing directory {0} does not seem to exist.".format(procdir)
        return 1
    slcdir = os.path.join(procdir,'SLC')
    if not os.path.exists:
        print >>sys.stderr, "\nERROR:"
        print >>sys.stderr, "SLC directory {0} does not seem to exist.".format(slcdir)
        return 1

############################################################ Check that existing directories are present
    masterslcdir = os.path.join(slcdir,masterdate.strftime('%Y%m%d'))
    if not os.path.exists:
        print >>sys.stderr, "\nERROR:"
        print >>sys.stderr, "Master SLC directory {0} does not seem to exist.".format(masterslcdir)
        return 1
    demdir = os.path.join(procdir,'DEM')
    if not os.path.exists:
        print >>sys.stderr, "\ERROR:"
        print >>sys.stderr, "DEM directory {0} does not seem to exist.".format(masterslcdir)
        return 1


############################################################ Parse multi look parameters from parameter files
    gc.rglks = int(grep('range_looks',os.path.join(masterslcdir,masterdate.strftime('%Y%m%d')+'.slc.mli.par')).split(':')[1].strip())
    gc.azlks = int(grep('azimuth_looks',os.path.join(masterslcdir,masterdate.strftime('%Y%m%d')+'.slc.mli.par')).split(':')[1].strip())
    # Set the output DEM resolution based on the number of looks of the existing data
    # Else force the old default value of the oversampling factor to be 2 as it was previously hardcoded.
    if gc.rglks == 20 and gc.azlks == 4:
        gc.outres = 0.001
    else:
        gc.outres = 0.0001
                
############################################################ Create a report file
    if not reportfile:
        reportfile = os.path.join(procdir,'{0}-coreg-report.txt'.format(framename))
    starttime = dt.datetime.now()
    with open(reportfile,'w') as f:
        f.write('LiCSAR_02_coreg processing started {0}.\n'.format(dt.datetime.now().strftime('%Y-%m-%d %H:%M')))
        f.write('\nMaster date: {0}\n'.format(masterdate))
    
############################################################ Geocode the DEM to master image
        if not ignrMstr:
        ############################################################Create Geoencoded directory
            geodir = os.path.join(procdir,'geo')
            if os.path.exists(geodir):
                shutil.rmtree(geodir)
            os.mkdir(geodir)
            rc = geocode_dem(masterslcdir,geodir,demdir,procdir,masterdate,gc.outres)
            with open(reportfile,'a') as f:
                if rc == 0:
                    f.write('\nMaster geocoding completed successfully\n')
                if rc == 1:
                    f.write('\nMaster geocoding encountered a problem during the lookup table creation')
                    return 1
                if rc == 2:
                    f.write('\nMaster geocoding encountered a problem during the amplitude simulation')
                    return 1
                if rc == 3:
                    f.write('\nMaster geocoding encountered a problem during the cross-correlation offset estimation')
                    return 1
                if rc == 4:
                    f.write('\nMaster geocoding encountered a problem during the refinement of the offsets')
                    return 1
                if rc == 5:
                    f.write('\nMaster geocoding encountered a problem during the creation of the DEM MLI')
                    return 1
                if rc == 6:
                    f.write('\nMaster geocoding encountered a problem during the creation of the master height file')
                    return 1
            
    ############################################################ Link master slc dir in resample dir
            rslcdir = os.path.join(procdir,'RSLC')
            if os.path.exists(rslcdir):
                shutil.rmtree(rslcdir)
            rc = link_master_rslc(masterslcdir,rslcdir,masterdate, lq, job_id)
            if rc > 0:
                f.write('\nProblem creating a link to master SLC directory in RSLC directory.')
                return 1
        else:
            f.write('\nMaster geocode ignore flag given - skipping geocoding')
            rslcdir = os.path.join(procdir,'RSLC')

    
############################################################ Create list of slave images
    if coregListFile:
        slavedatelist = list()
        with open(coregListFile,'r') as f:
            for line in f:
                slavedatelist.append(dt.datetime.strptime(line.strip(),'%Y%m%d'))

    else:
        slavelist = [s for s in os.listdir(slcdir) if s != masterdate.strftime('%Y%m%d')]
        slavedatelist = [dt.datetime.strptime(sd,'%Y%m%d') for sd in slavelist]

    timedelta = [sd.date()-masterdate for sd in slavedatelist] # difference between slaves and master dates
    sortix = np.argsort(np.abs(timedelta))
    slavedatelistsort = [slavedatelist[i] for i in sortix] #sorted differences
    
############################################################ Coregister slave images to master
    createdRslcs = []
    for sd in slavedatelistsort:
        lut = os.path.join(rslcdir,sd.strftime('%Y%m%d'),masterdate.strftime('%Y%m%d_')+sd.strftime('%Y%m%d.slc.mli.lt'))
        if os.path.exists(lut) and not forceRecreate:
            rc = recoreg_slave(sd,slcdir,rslcdir,masterdate,framename,procdir,lq)
            with open(reportfile,'a') as f:
                if rc == 0:
                    f.write('\nAcquisition {0} has been coregistered correctly.'.format(sd))
                    createdRslcs.append(sd)
                    if not keepSLCs:
                        print 'Removing SLC directory...'
                        imdir = os.path.join(slcdir,sd.strftime('%Y%m%d'))
                        shutil.rmtree(imdir)
                    if removeMosaic:
                        rslc = os.path.join(rslcdir,sd.strftime('%Y%m%d'),sd.strftime('%Y%m%d.rslc'))
                        os.remove(rslc)
                        f.write('\nRemoved mosaiced rslc for date {0}'.format(sd))
                elif rc == 1:
                    f.write('\nAcquisition {0} had a problem during the tab file creation.'.format(sd))
                elif rc == 3:
                    f.write('\nAcquisition {0} had a problem during the resampling'.format(sd))
                elif rc == 5:
                    f.write('\nAcquisition {0} had a problem during the recropping of the master'.format(sd))
                elif rc == 6: #These should not happen...
                    f.write('\nAcquisition {0} has not been coregistered'.format(sd))
                elif rc == 7:
                    f.write('\nAcquisition {0} does not have a lookup table'.format(sd))
        else:
            rc = coreg_slave(sd,slcdir,rslcdir,masterdate,framename,procdir, lq, job_id)
            with open(reportfile,'a') as f:
                if rc == 0:
                    f.write('\nAcquisition {0} has been coregistered correctly.'.format(sd))
                    createdRslcs.append(sd)
                    if not keepSLCs:
                        print 'Removing SLC directory...'
                        imdir = os.path.join(slcdir,sd.strftime('%Y%m%d'))
                        shutil.rmtree(imdir)
                    if removeMosaic:
                        rslc = os.path.join(rslcdir,sd.strftime('%Y%m%d'),sd.strftime('%Y%m%d.rslc'))
                        os.remove(rslc)
                        f.write('\nRemoved mosaiced rslc for date {0}'.format(sd))
                elif rc == 1:
                    f.write('\nAcquisition {0} had a problem during the tab file creation.'.format(sd))
                elif rc == 2:
                    f.write('\nAcquisition {0} had a problem during the geometric coregistration.'.format(sd))
                elif rc == 3:
                    f.write('\nAcquisition {0} had a problem during the resampling'.format(sd))
                elif rc == 4:
                    f.write('\nAcquisition {0} had a problem during the spectral diversity'.format(sd))
                elif rc == 5:
                    f.write('\nAcquisition {0} had a problem during the recropping of the master'.format(sd))
        
    
################################################################################
# Execute main function
################################################################################
if __name__ == "__main__":
    sys.exit(main())
