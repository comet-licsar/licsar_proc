#!/usr/bin/env python
"""
Sets up master image for further processing, without taking into account any slave images

========
Overview
========


=========
Changelog
=========
February 2019: Added option to automatically choose master (automaster) (Milan Lazecky, Uni of Leeds)
October 2017: Added option to include output resolution in decimal degrees (Emma Hatton, Uni of Leeds)
May 2017: Added step to geocode and tiff the inc and psi outputs (Emma Hatton, Uni of Leeds)
February 2017: Original implementation (Karsten Spaans, Uni of Leeds)


=====
Usage
=====

LiCSAR_setup_master -h 
Display this message

LiCSAR_setup_master.py -v
Display the version

LiCSAR_setup_master.py -f <frame_name> -d <output directory> [-m <master date yyyymmdd> -r <range looks> -a <azimuth looks> -j <job_id> -D <path to custom DEM (tif)> -A]
  (parameter -A or --automaster : will attempt to get a master date automatically)
  (parameter -V XXX would perform a dryrun, i.e. only check for possible appropriate masters, from the dataset since last XXX days)
  (parameter -e would mean allowing also latest images to become master - useful for earthquake responder)
  (parameter -D for custom DEM in tif format - e.g. download ALOS World 3D tiles and merge them by gdal_merge.py -o out.tif -of GTiff -a_nodata -32768 *.tif)
"""
################################################################################
# Imports
################################################################################
import getopt
import os
import shutil
import sys
import datetime as dt
from datetime import datetime, timedelta
import subprocess as subp
#Import global configuration, LiCSquery and gamma functions
import global_config as gc
import LiCSquery as lq
from gamma_functions import *
#Import functions from the LiCSAR processing chain
from LiCSAR_lib.mk_imag_lib import check_master_bursts, check_bursts, make_frame_image
from LiCSAR_lib.coreg_lib import link_master_rslc, geocode_dem
from LiCSAR_lib.LiCSAR_misc import Usage,cd
from LiCSAR_lib import s1data
import framecare as fc
import pandas as pd

################################################################################
# Master job started function
################################################################################
def master_job_started(job_id, acq_date):
    if acq_date and job_id != -1:
        acq_date = acq_date.strftime('%Y%m%d')
        # update jobs table with status and time
        lq.set_job_started(job_id)
        # update job requests table with status
        lq.set_job_request_started_master(job_id, acq_date)


################################################################################
# Master job finished clean function
################################################################################
def master_job_finished_clean(job_id, acq_date):
    if acq_date and job_id != -1:
        acq_date = acq_date.strftime('%Y%m%d')
        # Don't update the finished time and status of the job table, in case multiple dates being
        # processed, this will be done later by a separate python script called at the end of a batch
        # Update the job requests table to state the master finished cleanly
        lq.set_job_request_finished_clean(job_id, acq_date)


################################################################################
# master job finished failed function
################################################################################
def master_job_finished_failed(job_id, acq_date, status=-1):
    if acq_date and job_id != -1:
        acq_date = acq_date.strftime('%Y%m%d')
        # Don't update the finished time and status of the job table, in case multiple dates being
        # processed, this will be done later by a separate python script called at the end of a batch
        if status == -1:
            # Update the job requests table to state the master finished cleanly
            lq.set_job_request_finished_master_fail(job_id, acq_date)
        else:
            # Update the job requests table to state the master finished cleanly
            lq.set_job_request_finished_master_fail(job_id, acq_date, status)


################################################################################
# main function
################################################################################
def main(argv=None):
############################################################ Initialise parameters
    if argv == None:
        argv = sys.argv
    # Parameter initialisation and checking
    framename=[]
    burstidfile = []
    polygonfile = []
    procdir=[]
    startdate=[]
    enddate=[]
    masterdate=[]
    job_id = -1
    automaster = 0
    days_limit = 150
    days_limit_POD = 22
    customdem = ''
    dryrun = False
    autodownload = False
    
############################################################ Parse argument list
    try:
        try:
            opts, args = getopt.getopt(argv[1:], "vhf:V:d:j:m:a:r:o:A:D:e:", ["version", "help","automaster"])
        except getopt.error as msg:
            raise Usage(msg)
        for p, a in opts:
            if p == '-h' or p == '--help':
                print(__doc__)
                return 0
            elif p == '-v' or p == '--version':
                print("")
                print("Current version: %s" % gc.config['VERSION'])
                print("")
                return 0
            elif p == '-f':
                framename = a
            elif p == '-d':
                procdir = a
            elif p == '-D':
                customdem = a
            elif p == '-V':
                days_limit = int(a)
                dryrun = True
            elif p == '-j':
                job_id = int(a)
            elif p == '-m':
                masterdate = dt.date(int(a[:4]),int(a[4:6]),int(a[6:8]))
                autodownload = True
            elif p == '-A' or p == '--automaster':
                automaster = 1
            elif p == '-a':
                gc.azlks = int(a)
            elif p == '-r':
                gc.rglks = int(a)
            elif p == '-o':
                gc.outres = float(a)
            elif p == '-e':
                days_limit_POD = 0

        if not (framename or polygonfile or burstidfile):
            raise Usage('No frame given, please define the -f option!')

        if not procdir:
            raise Usage('No output data directory given, -d is not optional!')

        if job_id == -1:
            print("This processing is not outputting any products to the database.")

    except Usage as err:
        print("\nERROR:", file=sys.stderr)
        print("  "+str(err.msg), file=sys.stderr)
        print("\nFor help, use -h or --help.\n", file=sys.stderr)
        return 2

############################################################ Ensure connection to database
    if not lq.connection_established():
        print("\nERROR:", file=sys.stderr)
        print("Could not establish a stable database connection. No processing can happen.", file=sys.stderr)

        return 1

############################################################ Make sure processing directory exists
    if not os.path.exists(procdir):
        if not dryrun:
            os.mkdir(procdir)

############################################################ Check frame exists
    if framename:
        testlist = lq.check_frame(framename)
        if not testlist:
            print("\nERROR:", file=sys.stderr)
            print("Could not find frame {0} in database".format(framename), file=sys.stderr)
            return 1


############################################################ Get master and ensure no missing bursts
    # Check if master is in date list, and has all bursts. Missing bursts ok
    # for slaves, but not for master.
    if masterdate:
        #the timedelta below is a fix for midnight error in case of master:
        print('rechecking S1 data for this date')
        todown = s1data.check_and_import_to_licsinfo(framename,masterdate - timedelta(days=1), masterdate + timedelta(days=1), reingest = False)
        filelist = lq.get_frame_files_period(framename,masterdate-dt.timedelta(days=1),masterdate+dt.timedelta(days=1))
        midnighterror = False
        for f in filelist:
            if f[1] > masterdate: midnighterror = True
        burstlist = lq.get_bursts_in_frame(framename)
        rc = check_master_bursts(framename,burstlist,masterdate,[masterdate],lq, midnighterror)
        if rc != 0:
            return 1
    else:
        print('checking for S1 data not ingested to licsinfo db')
        todown = s1data.check_and_import_to_licsinfo(framename,(dt.date.today() - timedelta(days=days_limit)), reingest = False)
        print(todown)
        try:
            burstlist, filelist, dates = check_bursts(framename,dt.date(2015,10,0o1),dt.date.today(),lq)
        except:
            burstlist = None
            filelist = None
            dates = None
        if automaster and dates:
            print('doing automatic selection of master')
            #rc will be changed to 0 if a proper master is found and checked
            rc = 1
            frameburstlist = lq.get_bursts_in_frame(framename)
            framebursts = pd.DataFrame.from_records(frameburstlist)[0].to_list()
            cands = []
            for m in sorted(list(dates)):
                #master should be from files with POD (<3 weeks) and available (>90 days... or just 'days_limit')
                #but if we focus on earthquake response, we may use the latest ones also for master - so lets try
                if m > (dt.date.today() - timedelta(days=days_limit)) and m < (dt.date.today() - timedelta(days=days_limit_POD)) and rc == 1:
                    a = m.strftime('%Y%m%d')
                    masterdate = dt.date(int(a[:4]),int(a[4:6]),int(a[6:8]))
                    print('Checking {0} date as master'.format(masterdate))
                    filelist = lq.get_frame_files_period(framename,masterdate,masterdate)
                    rc = check_master_bursts(framename,frameburstlist,masterdate,[masterdate],lq)
                    if dryrun == True and rc == 0:
                        cands.append(str(masterdate))
                        #print('master candidates: '+str(masterdate))
                        rc = 1
                    if rc == 1 and not dryrun:
                        masterfiles = s1data.get_images_for_frame(framename, masterdate-timedelta(days=1), masterdate+timedelta(days=1))
                        if len(masterfiles) > len(filelist):
                            print('checking for possible cross-relorb problem')
                            filelist_pd = pd.DataFrame.from_records(filelist)
                            for mfile in masterfiles:
                                if mfile not in filelist_pd[2].to_list():
                                    print('new file {} identified using ASF - checking it'.format(mfile))
                                    mfile_bursts = lq.sqlout2list(lq.get_bursts_in_file(mfile))
                                    for mburst in mfile_bursts:
                                        if not mburst in framebursts:
                                            #print('checking burst {}'.format(mburst))
                                            if fc.check_and_fix_burst(mburst, framebursts):
                                                print('WARNING - burst definitions changed. if the init script fails, please rerun it')
            if dryrun:
                print('master candidates: ')
                print(cands)
                return 1
            if rc == 0:
                print('\nContinuing with the selected date {0} as master'.format(masterdate))
            else:
                if not dryrun:
                    print('\nAutomatic selection of master failed')
                return 1
        else:
            print('\nNo master date given. Please use the -m option to define one of these choices for the master:\n{0}'.format(', '.join([m.strftime('%Y%m%d') for m in sorted(list(dates)) if m != masterdate])), file=sys.stderr)
            return 1

############################################################ Update job Start
    # Log to DB that the processing has started
    master_job_started(job_id, masterdate)

############################################################ Get Frame image
    print('\nUnzipping, converting and merging image files for master {0}...'.format(masterdate))
    imburstlist = lq.get_frame_bursts_on_date(framename,masterdate)
    rc = make_frame_image(masterdate,framename,imburstlist,procdir, lq,job_id, autodownload = autodownload)
    if rc != 0:
        print('error initialising frame reference')
        return 1
    
    reportfile = os.path.join(procdir,'{0}-setup_master-report.txt'.format(framename))
    starttime = dt.datetime.now()

    # ELH - addition - output a file with polygon data
    polyfile = os.path.join(procdir,'{0}-poly.txt'.format(framename))
    with open(polyfile, 'w') as f:
        frame_poly = lq.get_polygon(framename)[0]
        frame_poly_cor = []
        for framepoly in frame_poly:
            if framepoly != None:
                frame_poly_cor.append(framepoly)
        frame_poly = frame_poly_cor
        frame_poly_zip = list(zip(frame_poly[::2], frame_poly[1::2]))
        for i in frame_poly_zip:
            f.write('{0} {1}\n'.format(i[0], i[1]))

############################################################ Create Report
    with open(reportfile,'w') as f:
        f.write('LiCSAR_setup_master processing started {0}.\n'.format(dt.datetime.now().strftime('%Y-%m-%d %H:%M')))
        f.write('\nMaster date: {0}\n'.format(masterdate))


    with open(reportfile,'a') as f:
        if rc == 0:
            f.write('Acquisition {0} completed the image setup successfully.'.format(masterdate))
        elif rc == 1:
            f.write('\nAcquisition {0} had a problem during the file reading.'.format(masterdate))
#                shutil.rmtree(datedir)
            master_job_finished_failed(job_id, masterdate)
            return 1
        elif rc == 2:
            f.write('\nAcquisition {0} had a problem during the copying and merging of bursts.'.format(masterdate))
#                shutil.rmtree(datedir)
            master_job_finished_failed(job_id, masterdate)
            return 1
        elif rc == 3:
            f.write('\nAcquisition {0} had a problem during the mosaicing of subswaths.'.format(masterdate))
#                shutil.rmtree(datedir)
            master_job_finished_failed(job_id, masterdate)
            return 1
        elif rc == 4:
            f.write('\nAcquisition {0} has missing bursts in critical location.'.format(masterdate))
#                if os.path.exists(datedir):
#                    shutil.rmtree(datedir)
            master_job_finished_failed(job_id, masterdate)
            return 1

############################################################ Crop the DEM
    print('\n\n')
    with cd(procdir):
        if customdem:
            demcall = 'LiCSAR_01_mk_crop_extDEM DEM/dem_crop 0 {}'.format(customdem)
        else:
            demcall = 'LiCSAR_01_mk_crop_extDEM DEM/dem_crop'
        os.system(demcall)

############################################################ Create Directory structure
    slcdir = os.path.join(procdir,'SLC')
    geodir = os.path.join(procdir,'geo')
    logdir = os.path.join(procdir,'log')
    if not os.path.exists(geodir):
        os.mkdir(geodir)
    if not os.path.exists(logdir):
        os.mkdir(logdir)
    demdir = os.path.join(procdir,'DEM')
    masterslcdir = os.path.join(slcdir,masterdate.strftime('%Y%m%d'))
    rslcdir = os.path.join(procdir,'RSLC')
    
############################################################ Geocode the master
    print('\n\n')
    rc = geocode_dem(masterslcdir,geodir,demdir,procdir,masterdate,gc.outres, skip_fit = False)
    if rc != 0:
        print('some error geocoding DEM. trying to skip the fitting (S1 orbits are very precise now)')
        rc = geocode_dem(masterslcdir,geodir,demdir,procdir,masterdate,gc.outres, skip_fit = True)
    with open(reportfile,'a') as f:
        if rc == 0:
            f.write('\nMaster geocoding completed successfully\n')
            gtcall = ['create_geoctiff_lookangles.sh', procdir, masterdate.strftime('%Y%m%d')]
            print('now do :')
            print(gtcall)
            #try:
            #    gt_code = subp.check_call(gtcall)
            #    if gt_code != 0:
            #        print('Something went wrong during the geotiff creation - non zero return.') 
            #except:
            #    print('Something went wrong during the geotiff creation - except call.')
            #    #f.write('\nSomething went wrong during the geotiff creation\n')
            #the line below should be removed. but i keep it, otherwise i would cause an error further in licsar_initiate_new_frame.sh
            #if gt_code != 0:
            #    print('Something went wrong during the geotiff creation - non zero return.') 
        if rc == 1:
            f.write('\nMaster geocoding encountered a problem during the lookup table creation')
            master_job_finished_failed(job_id, masterdate)
            return 1
        if rc == 2:
            f.write('\nMaster geocoding encountered a problem during the amplitude simulation')
            master_job_finished_failed(job_id, masterdate)
            return 1
        if rc == 3:
            f.write('\nMaster geocoding encountered a problem during the cross-correlation offset estimation')
            master_job_finished_failed(job_id, masterdate)
            return 1
        if rc == 4:
            f.write('\nMaster geocoding encountered a problem during the refinement of the offsets')
            master_job_finished_failed(job_id, masterdate)
            return 1
        if rc == 5:
            f.write('\nMaster geocoding encountered a problem during the creation of the DEM MLI')
            master_job_finished_failed(job_id, masterdate)
            return 1
        if rc == 6:
            f.write('\nMaster geocoding encountered a problem during the creation of the master height file')
            master_job_finished_failed(job_id, masterdate)
            return 1
############################################################ Update the resampled slc
    rc = link_master_rslc(masterslcdir,rslcdir,masterdate,lq, job_id)
    with open(reportfile,'a') as f:
        if rc > 0:
            f.write('\nProblem creating a link to master SLC directory in RSLC directory.')
            master_job_finished_failed(job_id, masterdate)
            return 1
############################################################ Update frames.csv file
    print('updating frames.csv file')
    gpan = fc.frame2geopandas(framename)
    fc.export_frames_to_licsar_csv(gpan, store_zero = True)
############################################################ fix potential cross-bursts
    print('performing check/fix of potential cross-track burst definitions')
    fc.check_and_fix_all_bursts_in_frame(framename)
############################################################ Update job database
    # Write to DB output the time the processing has finished and update its status.
    master_job_finished_clean(job_id, masterdate)
    return 0

################################################################################
# Execution main function
################################################################################
if __name__ == "__main__":
    sys.exit(main())
