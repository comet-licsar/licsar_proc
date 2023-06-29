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
2020+: some minor changes to reflect LiCSAR development (Milan Lazecky, Uni of Leeds)
August 2016: Original implementation (Karsten Spaans, Uni of Leeds)


=====
Usage
=====
LiCSAR_02_coreg.py -f <framename> -d </path/to/processing/location> -m <masterdate>

    -f    Name of the frame to be processed in database
    -d    Path to the processing location. 
    -m    Date of master image. If not given, it will give a list of available dates
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
    (-E   little tweaks for eq responder - e.g. skip SD quality check)

guys.. if you work with existing LiCSAR frame, please use -i to prevent re-geocoding

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
from gamma_functions import *

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
    cleandone = True
    job_id = -1
    coregListFile = None
    keepSLCs = False
    removeMosaic = False
    removeSSRslc = False
    reportfile = None
    forceRecreate = False
    eidp = False
############################################################ Parse arguments
    try:
        try:
            opts, args = getopt.getopt(argv[1:], "vhkiCcREl:f:d:m:j:y:z:p:T:", ["version", "help"])
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
            elif o == '-i':
                ignrMstr = True
            elif o == '-E':
                eidp = True
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
            print("This processing is not outputting any products to the database.")
############################################################ Setup processing database
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

############################################################ Ensure we have a connection to the database
    #Check if a DB connection can be established
    if not lq.connection_established():
        print("\nERROR:", file=sys.stderr)
        print("Could not establish a stable database connection. No processing can happen.", file=sys.stderr)

        return 1

############################################################ Check that parameters are in database
    if not lq.check_frame(framename):
        print("\nERROR:", file=sys.stderr)
        print("Frame {0} was not found in database.".format(framename), file=sys.stderr)

    if not os.path.exists(procdir):
        print("\nERROR:", file=sys.stderr)
        print("Processing directory {0} does not seem to exist.".format(procdir), file=sys.stderr)
        return 1
    slcdir = os.path.join(procdir,'SLC')
    if not os.path.exists(slcdir):
        print("\nERROR:", file=sys.stderr)
        print("SLC directory {0} does not seem to exist.".format(slcdir), file=sys.stderr)
        return 1

############################################################ Check that existing directories are present
    masterslcdir = os.path.join(slcdir,masterdate.strftime('%Y%m%d'))
    if not os.path.exists(masterslcdir):
        print("\nERROR:", file=sys.stderr)
        print("Master SLC directory {0} does not seem to exist.".format(masterslcdir), file=sys.stderr)
        return 1
    #demdir = os.path.join(procdir,'DEM')
    #if not os.path.exists:
    #    print("\ERROR:", file=sys.stderr)
    #    print("DEM directory {0} does not seem to exist.".format(masterslcdir), file=sys.stderr)
    #    return 1


############################################################ Parse multi look parameters from parameter files
    gc.rglks = int(grep1('range_looks',os.path.join(masterslcdir,masterdate.strftime('%Y%m%d')+'.slc.mli.par')).split(':')[1].strip())
    gc.azlks = int(grep1('azimuth_looks',os.path.join(masterslcdir,masterdate.strftime('%Y%m%d')+'.slc.mli.par')).split(':')[1].strip())
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
            demdir = os.path.join(procdir,'DEM')
            if not os.path.exists(demdir):
                print("\ERROR:", file=sys.stderr)
                print("DEM directory {0} does not seem to exist.".format(demdir), file=sys.stderr)
                return 1
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
        else:
            f.write('\nMaster geocode ignore flag given - skipping geocoding')

        ############################################################ Link master slc dir in resample dir, regenerate mosaic
        masterslcdir = os.path.join(procdir, 'SLC', masterdate.strftime('%Y%m%d'))
        mastermosaic = os.path.join(masterslcdir,masterdate.strftime('%Y%m%d')+'.slc')
        rslcdir = os.path.join(procdir,'RSLC')
        try:
            rc = link_master_rslc(masterslcdir, rslcdir, masterdate, lq, job_id=-1, overwrite = False)
        except:
            f.write('\nDEBUG: ERROR creating a link to master SLC directory in RSLC directory, but continuing.')
        if not os.path.exists(mastermosaic):
            print('regenerating master mosaic (some 2 minutes)')
            slctab = os.path.join(masterslcdir,masterdate.strftime('%Y%m%d')+'.mosaic.tab')
            logmosaic = os.path.join(masterslcdir,masterdate.strftime('%Y%m%d')+'.mosaic.log')
            iwfiles = glob(os.path.join(masterslcdir,masterdate.strftime('%Y%m%d')+'.IW*.slc'))
            swathlist = []
            for iwfile in iwfiles:
                swathlist.append(iwfile.split('/')[-1].split('.')[1])
            swathlist.sort()
            rc, msg = make_SLC_tab(slctab,mastermosaic,swathlist)
            rc = SLC_mosaic_S1_TOPS(slctab,mastermosaic,gc.rglks,gc.azlks,logmosaic)
            os.remove(slctab)
        rslcdir = os.path.join(procdir,'RSLC')
        if not os.path.exists(os.path.join(rslcdir, masterdate.strftime('%Y%m%d'), masterdate.strftime('%Y%m%d')+'.rslc.mli.par')):
            #shutil.rmtree(os.path.join(rslcdir, masterdate.strftime('%Y%m%d')))
            rc = link_master_rslc(masterslcdir,rslcdir,masterdate, lq, job_id)
            if rc > 0:
                f.write('\nProblem creating a link to master SLC directory in RSLC directory.')
                return 1
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
    lutdir = os.path.join(procdir,'LUT')  # changed to expected place
    for sd in slavedatelistsort:
        lut1 = os.path.join(rslcdir,sd.strftime('%Y%m%d'),masterdate.strftime('%Y%m%d_')+sd.strftime('%Y%m%d.slc.mli.lt'))
        lut2 = os.path.join(lutdir,sd.strftime('%Y%m%d'),masterdate.strftime('%Y%m%d_')+sd.strftime('%Y%m%d.slc.mli.lt'))
        if (os.path.exists(lut1) or os.path.exists(lut2)) and not forceRecreate:
            rc = recoreg_slave(sd,slcdir,rslcdir,masterdate,framename,procdir,lq)
            with open(reportfile,'a') as f:
                if rc == 0:
                    f.write('\nAcquisition {0} has been coregistered correctly.'.format(sd))
                    createdRslcs.append(sd)
                    if not keepSLCs:
                        print('Removing SLC directory...')
                        imdir = os.path.join(slcdir,sd.strftime('%Y%m%d'))
                        if os.path.exists(imdir):
                            shutil.rmtree(imdir)
                    if removeMosaic:
                        rslc = os.path.join(rslcdir,sd.strftime('%Y%m%d'),sd.strftime('%Y%m%d.rslc'))
                        if os.path.exists(rslc):
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
                    f.write('\nAcquisition {0} does not have a lookup table or other LUT-related error'.format(sd))
        else:
            rc = coreg_slave(sd,slcdir,rslcdir,masterdate,framename,procdir, lq, job_id, eidp = eidp)
            with open(reportfile,'a') as f:
                if rc == 0:
                    f.write('\nAcquisition {0} has been coregistered correctly.'.format(sd))
                    createdRslcs.append(sd)
                    if not keepSLCs:
                        print('Removing SLC directory...')
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
        #if rc == 0:
        if cleandone:
            slavestr = sd.strftime('%Y%m%d')
            masterstr = masterdate.strftime('%Y%m%d')
            slaverslcdir = os.path.join(rslcdir,sd.strftime('%Y%m%d'))
            #for ext in ['']
            os.system('rm {0}/{1}.IW?.rslc.[0-9]*.[0-9]* 2>/dev/null'.format(slaverslcdir,slavestr))
            os.system('rm {0}_{1}.* {1}.* SLC_interp_lt_ScanSAR.* core.* 2>/dev/null'.format(masterstr,slavestr))
    try:
        lq.close_db_and_tunnel()
    except:
        print('')
################################################################################
# Execute main function
################################################################################
if __name__ == "__main__":
    sys.exit(main())
