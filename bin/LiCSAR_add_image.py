#!/usr/bin/env python
"""
Form interferogram for newly acquired images
"""

################################################################################
#Imports
################################################################################
import getopt
import os
import sys
import subprocess as subp
import numpy as np
import datetime as dt
import shutil
# import global config, LiCSquer and gamme functions
import global_config as gc
import LiCSquery as lq
from gamma_functions import *
#Import specfic LiCSAR functions
from LiCSAR_lib.mk_imag_lib import make_frame_image, check_missing_bursts
from LiCSAR_lib.coreg_lib import coreg_slave, link_master_rslc
from LiCSAR_lib.ifg_lib import make_interferogram
from LiCSAR_lib.unwrp_lib import do_unwrapping
from LiCSAR_lib.LiCSAR_misc import Usage

################################################################################
# Standard job started function
################################################################################
def standard_job_started(job_id, acq_date):
    if acq_date and job_id != -1:
        acq_date = acq_date.strftime('%Y%m%d')
        # update jobs table with status and time
        lq.set_job_started(job_id)
        # update job requests table with status
        lq.set_job_request_started_standard(job_id, acq_date)

################################################################################
# Standard job finished clean function
################################################################################
def standard_job_finished_clean(job_id, acq_date):
    if acq_date and job_id != -1:
        acq_date = acq_date.strftime('%Y%m%d')
        # Don't update the finished time and status of the job table, in case multiple dates being
        # processed, this will be done later by a separate python script called at the end of a batch
        # Update the job requests table to state the master finished cleanly
        lq.set_job_request_finished_clean(job_id, acq_date)

################################################################################
# standard job finished failed function
################################################################################
def standard_job_finished_failed(job_id, acq_date, status=-1):
    if acq_date and job_id != -1:
        acq_date = acq_date.strftime('%Y%m%d')
        # Don't update the finished time and status of the job table, in case multiple dates being
        # processed, this will be done later by a separate python script called at the end of a batch

        if status == -1:
            # Update the job requests table to state the master finished cleanly
            lq.set_job_request_finished_standard_fail(job_id, acq_date)
        else:
            # Update the job requests table to state the master finished cleanly
            lq.set_job_request_finished_standard_fail(job_id, acq_date, status)

################################################################################
# main function
################################################################################
def main(argv=None):
############################################################ Init parameters
    if argv == None:
        argv = sys.argv

    framename = []
    procdir = []
    slavedate = []
    job_id = -1

############################################################ Parse args
    try:
        try:
            opts, args = getopt.getopt(argv[1:], "vhf:d:s:j:r:a:", ["version", "help"])
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
            elif o == '-j':
                job_id = int(a)
            elif o == '-s':
                slavedate = dt.date(int(a[:4]),int(a[4:6]),int(a[6:8]))
            elif o == '-a':
                gc.azlks = int(a)
            elif o == '-r':
                gc.rglks = int(a)
        if not framename:
            raise Usage('No frame, polygon or set of burst ids given, please use either the -f')
        if not procdir:
            raise Usage('No output data directory given, -d is not optional!')
        if not slavedate:
            raise Usage('No image date given, -s is not optional!')
        if job_id == -1:
            print("This processing is not outputting any products to the database.")


    except Usage as err:
        print("\nERROR:", file=sys.stderr)
        print("  "+str(err.msg), file=sys.stderr)
        print("\nFor help, use -h or --help.\n", file=sys.stderr)

        # Log failed processing
        standard_job_finished_failed(job_id, slavedate)

        return 2

############################################################ Ensure databade connection
    if not lq.connection_established():
        print("\nERROR:", file=sys.stderr)
        print("Could not establish a stable database connection. No processing can happen.", file=sys.stderr)

        return 1

############################################################ Check processed files exist
    if not os.path.exists(procdir):
        print("\nERROR:", file=sys.stderr)
        print("Processing directory {0} does not seem to exist.".format(procdir), file=sys.stderr)

        # Log failed processing
        standard_job_finished_failed(job_id, slavedate)
        return 1

    #geocoded directory
    geodir = os.path.join(procdir,'geo')
    if not os.path.exists(geodir):
        print("\nERROR:", file=sys.stderr)
        print("Geocoding directory {0} does not seem to exist.".format(geodir), file=sys.stderr)

        # Log failed processing
        standard_job_finished_failed(job_id, slavedate)
        return 1

    #get master name from height file
    masterdatestr = []
    geodirlist = os.listdir(geodir)
    for l in geodirlist:
        if l[-4:] == '.hgt' and l[0] == '2':
            masterdatestr = l[:8]
            
    if masterdatestr:
        masterdate = dt.date(int(masterdatestr[:4]),int(masterdatestr[4:6]),int(masterdatestr[6:8]))
    else:
        print("\nERROR:", file=sys.stderr)
        print("Could not find master date in geocoding directory {0}.".format(geodir), file=sys.stderr)

        # Log failed processing
        standard_job_finished_failed(job_id, slavedate)

        return 1

    # Check RSLC directory
    masterslcdir = os.path.join(procdir,'SLC',masterdatestr)
    masterrslcdir = os.path.join(procdir,'RSLC',masterdatestr)
    rslcdir = os.path.join(procdir,'RSLC')
    if not os.path.exists(masterrslcdir):
        rc = link_master_rslc(masterslcdir,rslcdir,masterdate,lq, job_id)
        if rc > 0:
            f.write('\nProblem creating a link to master SLC directory in RSLC directory.')
            return 1

############################################################ Get frame name
    if framename:
        testlist = lq.check_frame(framename)
        if not testlist:
            print("\nERROR:", file=sys.stderr)
            print("Could not find frame {0} in database".format(framename), file=sys.stderr)

            # Log failed processing
            standard_job_finished_failed(job_id, slavedate)

            return 1
        else:
            burstlist = lq.get_bursts_in_frame(framename)

############################################################ Get slave date data
    res = lq.get_frame_files_date(framename,slavedate)
    if not res:
        print("\nERROR:", file=sys.stderr)
        print("Could not find image date {0} in frame {1}".format(slavedate,framename), file=sys.stderr)

        # Log failed processing
        standard_job_finished_failed(job_id, slavedate)

        return 1

############################################################ Set job started
    standard_job_started(job_id, slavedate)

############################################################ Make the frame image from burst data
    imburstlist = lq.get_frame_bursts_on_date(framename,slavedate)
    missingbursts = [b for b in burstlist if not b in imburstlist]
    if len(missingbursts) < 1:
        # All bursts there, no problem
        print("All bursts for frame {0} seem to be have been acquired on {1}...".format(framename,slavedate))
        rc = make_frame_image(slavedate,framename,imburstlist,procdir, lq,job_id)
    else:
        # Missing one or more bursts, checking if in problematic loc
        print("One of more  bursts for frame {0} have not been acquired on the date {1}. Missing bursts: {2}".format(framename,slavedate,''.join(['\n'+m[0] for m in missingbursts])))
        print("Checking where missing bursts are.")
        if check_missing_bursts(burstlist,missingbursts):
            # Missing bursts are at edges, not a problem
            print("Missing bursts are not in critical location, continuing processing image {0}.".format(slavedate))
            rc = make_frame_image(slavedate,framename,imburstlist,procdir, lq,job_id)
        else:
            # Missing bursts in the middle, we're stuffed for 
            # this date, skip to next
            print("Missing bursts are in critical location, continuing processing with next image...")

            # Log failed processing
            standard_job_finished_failed(job_id, slavedate)

            return 1
    datedir = os.path.join(procdir,'SLC',slavedate.strftime('%Y%m%d'))
    if rc == 1:
        if os.path.exists(datedir):
            shutil.rmtree(datedir)

        # Log failed processing
        standard_job_finished_failed(job_id, slavedate)

        return 1
    elif rc == 2:
        shutil.rmtree(datedir)

        # Log failed processing
        standard_job_finished_failed(job_id, slavedate)

        return 1
    elif rc == 3:
        shutil.rmtree(datedir)

        # Log failed processing
        standard_job_finished_failed(job_id, slavedate)

        return 1
    elif rc == 4:
        if os.path.exists(datedir):
            shutil.rmtree(datedir)

        # Log failed processing
        standard_job_finished_failed(job_id, slavedate)

        return 1
############################################################ Coregester slave
    slcdir = os.path.join(procdir,'SLC')
    rslcdir = os.path.join(procdir,'RSLC')
    rc = coreg_slave(dt.datetime.combine(slavedate,dt.time()),slcdir,rslcdir,masterdate,framename,procdir, lq,job_id)
    if rc == 0:
        print('Removing SLC directory...')
        imdir = os.path.join(slcdir,slavedate.strftime('%Y%m%d'))
        shutil.rmtree(imdir)
    if rc == 1:
        imdir = os.path.join(rslcdir,slavedate.strftime('%Y%m%d'))
        if os.path.exists(imdir):
            shutil.rmtree(imdir)
        return 1
    if rc == 2:
        imdir = os.path.join(rslcdir,slavedate.strftime('%Y%m%d'))
        if os.path.exists(imdir):
            shutil.rmtree(imdir)
        return 1
    if rc == 3:
        imdir = os.path.join(rslcdir,slavedate.strftime('%Y%m%d'))
        if os.path.exists(imdir):
            shutil.rmtree(imdir)
        return 1
    if rc == 4:
        imdir = os.path.join(rslcdir,slavedate.strftime('%Y%m%d'))
        if os.path.exists(imdir):
            shutil.rmtree(imdir)
        return 1
    if rc == 5:
        imdir = os.path.join(rslcdir,slavedate.strftime('%Y%m%d'))
        if os.path.exists(imdir):
            shutil.rmtree(imdir)
        return 1
############################################################ Create unwrapped interferograms
    if len(missingbursts) > 1:
        cropgeodir = os.path.join(rslcdir,slavedate.strftime('%Y%m%d'),'geo')
        if os.path.exists(cropgeodir):
            shutil.rmtree(cropgeodir)
    procslavelist = sorted([d for d in os.listdir(rslcdir) if len(d) == 8 and d[0] == '2'])
    procslavelist_dt = [dt.date(int(d[:4]),int(d[4:6]),int(d[6:])) for d in procslavelist if dt.date(int(d[:4]),int(d[4:6]),int(d[6:])) != slavedate]

    ifg_fail = -1
    timediff = [pd-slavedate for pd in procslavelist_dt]
    #if new image is the latest date
    if timediff[-1] < dt.timedelta(0):
        # New image was acquired after all other acquisitions
        c = 1

        # count of how many IFGs failed
        ifg_fail = 0

        # While loop to ensure only the three closest combinations are formed
        while c < 4 and len(timediff) > 0:
            timediff.pop()
            masterthis = procslavelist_dt.pop()
            #Make the interferogram
            rc = make_interferogram(masterdate,masterthis,slavedate,procdir,lq, job_id)
            if rc == 1:
                print('\nAcquisition {0} had a problem during the offset file creation.'.format(slavedate))
                ifg_fail += 1
            if rc == 2:
                print('\nAcquisition {0} had a problem during the topographic phase estimation.'.format(slavedate))
                ifg_fail += 1
            if rc == 3:
                print('\nAcquisition {0} had a problem during the interferogram formation.'.format(slavedate))
                ifg_fail += 1
            if rc == 4:
                print('\nAcquisition {0} had a problem during the sunraster preview creation.'.format(slavedate))
                ifg_fail += 1
            if rc == 5:
                print('\nAcquisition {0} had a problem dhttp://www.leeds.ac.uk/uring the coherence estimation.'.format(slavedate))
                ifg_fail += 1

            if rc == 0:
                ifg = '{0}_{1}'.format(masterthis.strftime('%Y%m%d'), slavedate.strftime('%Y%m%d'))
                gtcall = ['create_geoctiff_single.sh', procdir, masterdatestr, ifg]
                print(gtcall)
                try:
                    gt_code = subp.check_call(gtcall)
                    if gt_code != 0:
                        print('Something went wrong during the geotiff creation - non zero return')
                except:
                    print('Something went wrong during the geotiff creation - except call.')


                # Interferogram successfully formed, start unwrapping
                ifgdir = os.path.join(procdir,'IFG')
                rc = do_unwrapping(masterdate.strftime('%Y%m%d'),ifg,ifgdir,procdir,lq,job_id)
                if rc ==1:
                    # Filtering
                    print('\Interferogram {0} had a problem during the filtering.'.format(ifg))
                elif rc == 2:
                    # Unwrapping
                    print('\Interferogram {0} had a problem during the unwrapping.'.format(ifg))
                if rc == 0 :
                    gtcall = ['create_geoctiff_unw.sh', procdir, masterdatestr, ifg]
                    print(gtcall)
                    try:
                        gt_code = subp.check_call(gtcall)
                    except:
                        print('Something went wrong during the unw geotiff creation - except call.')
                        if gt_code != 0:
                            print('Something went wrong during the unw geotiff creation - non zero return.')
                            
            ################################################
            # Deleting earliest slave                      #
            # Quick fix, needs to be done better probably! #
            ################################################
        


            c+=1

        if job_id != -1:
            rmdir = os.path.join(rslcdir,masterthis.strftime('%Y%m%d'))
            for part in ['.IW1.','.IW2.','.IW3.','.']:
                filethis = os.path.join(rmdir,'{0}{1}rslc'.format(masterthis.strftime('%Y%m%d'),
                                                                  part))
                os.remove(filethis)
                                   

    # is the date older than all others?
    elif timediff[0] > dt.timedelta(0):
        c = 1

        ifg_fail = 0

        # Go through and find the closest 3 dates and use this one as master
        while c < 4 and len(timediff) > 0:
            timediff.pop(0)
            slavethis = procslavelist_dt.pop(0)
            # New image (slavedate) is now master!
            rc = make_interferogram(masterdate,slavedate,slavethis,procdir,lq,job_id)
            if rc == 1:
                print('\nAcquisition {0} had a problem during the offset file creation.'.format(slavedate))
                ifg_fail += 1
            if rc == 2:
                print('\nAcquisition {0} had a problem during the topographic phase estimation.'.format(slavedate))
                ifg_fail += 1
            if rc == 3:
                print('\nAcquisition {0} had a problem during the interferogram formation.'.format(slavedate))
                ifg_fail += 1
            if rc == 4:
                print('\nAcquisition {0} had a problem during the sunraster preview creation.'.format(slavedate))
                ifg_fail += 1
            if rc == 5:
                print('\nAcquisition {0} had a problem during the coherence estimation.'.format(slavedate))
                ifg_fail += 1
            if rc == 0:
                ifg = '{0}_{1}'.format(masterthis.strftime('%Y%m%d'), slavedate.strftime('%Y%m%d'))
                gtcall = ['create_geoctiff_single.sh', procdir, masterdatestr, ifg]
                print(gtcall)
                try:
                    gt_code = subp.check_call(gtcall)
                    if gt_code != 0:
                        print('Something went wrong during the geotiff creation - nonzero return')
                except:
                    print('Something went wrong during the geotiff creation - except call.')

                # Interferogram successfully formed, start unwrapping
                ifgdir = os.path.join(procdir,'IFG')
                rc = do_unwrapping(masterdate.strftime('%Y%m%d'),ifg,ifgdir,procdir,lq,job_id)
                if rc ==1:
                    # Filtering
                    print('\Interferogram {0} had a problem during the filtering.'.format(ifg))
                elif rc == 2:
                    # Unwrapping
                    print('\Interferogram {0} had a problem during the unwrapping.'.format(ifg))
                if rc == 0 :
                    gtcall = ['create_geoctiff_unw.sh', procdir, masterdatestr, ifg]
                    print(gtcall)
                    try:
                        gt_code = subp.check_call(gtcall)
                    except:
                        print('Something went wrong during the unw geotiff creation - except call.')
                        if gt_code != 0:
                            print('Something went wrong during the unw geotiff creation - non zero return.')
                    
                    
            c+=1
        #####################
        # TODO Remove RSLCs #
        #####################
    else:
        # New image was acquired somewhere in the middle of the stack
        # Start with images acquired before current image
        timediff = np.array(timediff)
        ixafter = timediff < dt.timedelta(0)
        procslavelist_dt = np.array(procslavelist_dt)
        timediffafter = timediff[ixafter].tolist()
        procslavelistafter = procslavelist_dt[ixafter].tolist()
        c = 1

        # count of how many IFGs failed
        ifg_fail = 0

        # While loop to ensure only the three closest combinations are formed
        while c < 4 and len(timediffafter) > 0:
            timediffafter.pop()
            masterthis = procslavelistafter.pop()
            #make interferogram
            rc = make_interferogram(masterdate,masterthis,slavedate,procdir,lq,job_id)
            if rc == 1:
                print('\nAcquisition {0} had a problem during the offset file creation.'.format(slavedate))
                ifg_fail += 1
            if rc == 2:
                print('\nAcquisition {0} had a problem during the topographic phase estimation.'.format(slavedate))
                ifg_fail += 1
            if rc == 3:
                print('\nAcquisition {0} had a problem during the interferogram formation.'.format(slavedate))
                ifg_fail += 1
            if rc == 4:
                print('\nAcquisition {0} had a problem during the sunraster preview creation.'.format(slavedate))
                ifg_fail += 1
            if rc == 5:
                print('\nAcquisition {0} had a problem during the coherence estimation.'.format(slavedate))
                ifg_fail += 1
            if rc == 0:
                ifg = '{0}_{1}'.format(masterthis.strftime('%Y%m%d'), slavedate.strftime('%Y%m%d'))
                gtcall = ['create_geoctiff_single.sh', procdir, masterdatestr, ifg]
                print(gtcall)
                try:
                    gt_code = subp.check_call(gtcall)
                    if gt_code != 0:
                        print('Something went wrong during the geotiff creation - nonzero return.')
                except:
                    print('Something went wrong during the geotiff creation - except call.')


                # Interferogram successfully formed, start unwrapping
                ifgdir = os.path.join(procdir,'IFG')
                rc = do_unwrapping(masterdate.strftime('%Y%m%d'),ifg,ifgdir,procdir,lq,job_id)
                if rc ==1:
                    # Filtering
                    print('\Interferogram {0} had a problem during the filtering.'.format(ifg))
                elif rc == 2:
                        # Unwrapping
                    print('\Interferogram {0} had a problem during the unwrapping.'.format(ifg))
                if rc == 0 :
                    gtcall = ['create_geoctiff_unw.sh', procdir, masterdatestr, ifg]
                    print(gtcall)
                    try:
                        gt_code = subp.check_call(gtcall)
                    except:
                        print('Something went wrong during the unw geotiff creation - except call.')
                        if gt_code != 0:
                            print('Something went wrong during the unw geotiff creation - non zero return.')
            c+=1
        # Now images acquired after current image
        ixbefore = timediff > dt.timedelta(0)
        timediffbefore = timediff[ixbefore].tolist()
        procslavelistbefore = procslavelist_dt[ixbefore].tolist()
        # reset image counter, but not failed ifg counter
        c=1

        while c < 4 and len(timediffbefore) > 0:
            timediffbefore.pop(0)
            slavethis = procslavelistbefore.pop(0)
            # New image (slavedate) is now master!
            rc = make_interferogram(masterdate,slavedate,slavethis,procdir,job_id)
            if rc == 1:
                print('\nAcquisition {0} had a problem during the offset file creation.'.format(slavedate))
                ifg_fail += 1
            if rc == 2:
                print('\nAcquisition {0} had a problem during the topographic phase estimation.'.format(slavedate))
                ifg_fail += 1
            if rc == 3:
                print('\nAcquisition {0} had a problem during the interferogram formation.'.format(slavedate))
                ifg_fail += 1
            if rc == 4:
                print('\nAcquisition {0} had a problem during the sunraster preview creation.'.format(slavedate))
                ifg_fail += 1
            if rc == 5:
                print('\nAcquisition {0} had a problem during the coherence estimation.'.format(slavedate))
                ifg_fail += 1
            if rc == 0:
                ifg = '{0}_{1}'.format(masterthis.strftime('%Y%m%d'), slavedate.strftime('%Y%m%d'))
                gtcall = ['create_geoctiff_single.sh', procdir, masterdatestr, ifg]
                print(gtcall)
                try:
                    gt_code = subp.check_call(gtcall)
                    if gt_code != 0:
                        print('Something went wrong during the geotiff creation - nonzero return')
                except:
                    print('Something went wrong during the geotiff creation - except call.')
                    
                # Interferogram successfully formed, start unwrapping
                ifgdir = os.path.join(procdir,'IFG')
                rc = do_unwrapping(masterdate.strftime('%Y%m%d'),ifg,ifgdir,procdir,lq,job_id)
                if rc ==1:
                    # Filtering
                    print('\Interferogram {0} had a problem during the filtering.'.format(ifg))
                elif rc == 2:
                    # Unwrapping
                    print('\Interferogram {0} had a problem during the unwrapping.'.format(ifg))
                if rc == 0 :
                    gtcall = ['create_geoctiff_unw.sh', procdir, masterdatestr, ifg]
                    print(gtcall)
                    try:
                        gt_code = subp.check_call(gtcall)
                    except:
                        print('Something went wrong during the unw geotiff creation - except call.')
                        if gt_code != 0:
                            print('Something went wrong during the unw geotiff creation - non zero return.')
            c += 1
        #####################
        # TODO Remove RSLCs #
        #####################

    # ELH - Addition to ensure we have look angles for all files which are added to existing masters
    # Assume if one exists so does the other.
    psifile = ''.join([masterdatestr, '.geo.psi.tif'])
    if not os.path.join(procdir, 'GEOC/lookangles', psifile):
        gtcall = ['create_geoctiff_lookangles.sh', procdir, masterdatestr]
        try:
            gt_code = subp.check_call(gtcall)
            if gt_code != 0:
                print('Something went wrong during the look angle geotiff creation - non zero return.')                
        except:
            print('Something went wrong during the look angle geotiff creation - except call.')
    # Also make sure we have a polygon file
    polyfile = os.path.join(procdir,'{0}-poly.txt'.format(framename))
    if not polyfile:
        with open(polyfile, 'w') as f:
            frame_poly = lq.get_polygon(framename)[0]
            frame_poly_zip = list(zip(frame_poly[::2], frame_poly[1::2]))
            for i in frame_poly_zip:
                f.write('{0} {1}\n'.format(i[0], i[1]))

    # End of addition

############################################################ Update database on job issues
    if ifg_fail == -1:
        # Log failed processing, don't think it should be able to here!
        standard_job_finished_failed(job_id, slavedate)
    elif ifg_fail == 0:
        # Job finished cleanly output to DB
        standard_job_finished_clean(job_id, slavedate)
    elif ifg_fail > 0:
        # Log failed processing, some ifgs failed to be created
        # status code = 5000 (ifgs failed) + count of how many failed.
        status_code = 5000+ifg_fail
        standard_job_finished_failed(job_id, slavedate, status_code)
    else:
        # Log system failed, as this should be impossible
        print("\n\n###############################" \
              "\n\nERROR: LiCSAR_add_images.py has reached what should be an impossible condition " \
              "\n\nPlease contact the dev team to resolve.\n\n###############################\n\n")

    if ifg_fail != 0:
        return 1
    else:
        return 0

################################################################################
#  Execute main function
################################################################################
if __name__ == "__main__":
    sys.exit(main())
    
