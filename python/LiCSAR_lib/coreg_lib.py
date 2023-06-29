################################################################################
#Import Modules
################################################################################
import getopt
import os
import shutil
import re
import datetime
import sys
import numpy as np
import datetime as dt
import subprocess as subp
from LiCSAR_lib.LiCSAR_misc import *
import pandas as pd

from glob import glob

#Import Gamma functions
from gamma_functions import *
#And import the global configuration
import global_config as gc

################################################################################
#Link master rslc function
################################################################################
def link_master_rslc(masterslcdir, rslcdir, masterdate, lq, job_id, overwrite = True):
    """ Creates RSLC dir and links master SLC images
    """
############################################################Create directories
    rslcdir = os.path.realpath(rslcdir)
    masterslcdir = os.path.realpath(masterslcdir)
    if not os.path.exists(rslcdir):
        os.mkdir(rslcdir)
    masterrslcdir = os.path.join(rslcdir,masterdate.strftime('%Y%m%d'))
    if overwrite:
        if os.path.exists(masterrslcdir):
            shutil.rmtree(masterrslcdir)
    if not os.path.exists(masterrslcdir):
        os.mkdir(masterrslcdir)
############################################################ Iterate through 
                                                        #and syslink master slcs
    filelist = [x for x in os.listdir(masterslcdir) if '{0}.'.format(
                                    masterdate.strftime('%Y%m%d')) in x]
    filelistnew = [x.replace('.slc','.rslc') for x in filelist]
    try:
        for f1,f2 in zip(filelist,filelistnew):
            if not os.path.exists(os.path.join(masterrslcdir,f2)):
                os.symlink(os.path.join(masterslcdir,f1),
                           os.path.join(masterrslcdir,f2))
    except:
        print('Something went wrong linking the master SLCs to the '\
                                                    'RSLC directory.')
        return 1
############################################################ Update the job database
    # Populate the RSLC table with the location of the master RSLC
    if job_id != -1:
        print('Attempting to input new rslc MASTER product')
        lq.set_new_rslc_product(job_id, masterdate.strftime('%Y%m%d'),
                -1, masterdate.strftime('%Y%m%d')+'.rslc', masterrslcdir)
    return 0

################################################################################
#geocode dem function
################################################################################
def geocode_dem(masterslcdir,geodir,demdir,procdir,masterdate,outres = gc.outres, skip_fit = False):
    """ Geocodes the DEM to the master image geometry
    """
############################################################ Create radar coded DEM
    print("\nCreating radar coded DEM:")
    try:
        masterstr = masterdate.strftime('%Y%m%d')
    except:
        masterstr = masterdate
    logfilename = os.path.join(procdir,'log','gc_map_{0}.log'.format(masterstr))
    #create mli parameter file path, dem file path, dem parameter path
    mlipar = os.path.join(masterslcdir,masterstr+'.slc.mli.par')
    if not os.path.exists(mlipar):
        mlipar = os.path.join(masterslcdir,masterstr+'.rslc.mli.par')
    dem = os.path.join(demdir,'dem_crop.dem')
    dempar = os.path.join(demdir, 'dem_crop.dem_par')
    # Get dem res N and E (which are latitude and longitude resolutions?)
    with open(dempar, 'r') as f:
        for line in f:
            if 'post_lat' in line:
                demresN = float(line.split()[1])
            if 'post_lon' in line:
                demresE = float(line.split()[1])
    #calculate oversampling factors for lon/lat
    ovrfactN = str(-1.0*(demresN/outres))
    ovrfactE = str((demresE/outres))
    #create a dem segment file path? map segment?
    demseg = os.path.join(geodir,'EQA.dem')
    #look up table path?
    lut = os.path.join(geodir,masterstr+'.lt')
    #simulated SAR backscatter image
    simsar = os.path.join(geodir,masterstr+'.sim_sar')
    #zenith angle of surface normal vector n 
    u = os.path.join(geodir,'u')
    #orientation angle of n 
    v = os.path.join(geodir,'v')
    #local incidence angle
    inc = os.path.join(geodir,'inc')
    #projection angle
    psi = os.path.join(geodir,'psi')
    #pixel area normalization factor
    pix = os.path.join(geodir,'pix')
    #layover and shadow map
    lsmap = os.path.join(geodir,'ls_map')
############################################################ Calculate look up table
    print("Calculating look-up table...")
    #create a geocoded lookup table. Does this estimate demseg par?
    if not gc_map(mlipar,'-',dem,demseg,lut,ovrfactN,ovrfactE,simsar,
                                      u,v,inc,psi,pix,lsmap,'8','2','',
                                      logfilename):
        print("\nError:", file=sys.stderr)
        print("Something went wrong during the DEM lookup table creation.", file=sys.stderr)
        return 1 
############################################################ Simulate DEM amplitude
    mli2 = mlipar[:-4]
    if not skip_fit:
        #sigma 0 normalization area
        pixsigma = os.path.join(geodir,'pix_sigma0')
        #gamma 0 normalization area
        pixgamma = os.path.join(geodir,'pix_gamma0')
        logfilename = os.path.join(procdir,'log','pixel_area_{0}.log'.format(masterstr))
        print('Simulating DEM amplitude...')
        if not pixel_area(mlipar,demseg,lut,lsmap,inc,pixsigma,pixgamma,logfilename):
            # Error in amplitude simulation
            print("\nError:", file=sys.stderr)
            print("Something went wrong during the DEM amplitude simulation.", file=sys.stderr)
            return 2
    ############################################################ Create a Diff param file
    #    masterstr = masterdate.strftime('%Y%m%d')
        #the pixsigma file from the previous step
        mli1 = os.path.join(geodir,'pix_sigma0')
        #mli2 = mlipar[:-4] #The original master mli
        #Diff param file path
        diffpar = os.path.join(geodir,masterstr+'.diff_par')
        logfile = os.path.join(procdir,'log',
                               'create_diff_par_{0}.log'.format(masterstr))
        if not create_diff_par(mli2+'.par','-',diffpar,'1','0',logfile):
            print("\nError:", file=sys.stderr)
            print("Something went wrong during the DIFF parameter file creation", file=sys.stderr)
            return 3
    ############################################################ Estimate cross correlation
                                                    #offsets between sim and master original
        print('Estimating offsets...')
        #offset file
        offfile = os.path.join(geodir,masterstr+'.offs')
        #redundent variable?
        cofffile = os.path.join(geodir,'coffs')
        #patch cross-correlation
        ccpfile = os.path.join(geodir,masterstr+'.ccp')
        #text file containing offsets
        offsets = os.path.join(geodir,'offsets')
        logfile = os.path.join(procdir,'log',
                               'offset_pwrm_{0}.log'.format(masterstr))
        if not offset_pwrm(mli1,mli2,diffpar,offfile,ccpfile,'256','256',
                           offsets,'2','64','64','0.2',logfile):
            print("\nError:", file=sys.stderr)
            print("Something went wrong during the cross correlation"\
                                                                "offset estimation.", file=sys.stderr)
            return 3
    ############################################################ Fit offset function
        coffs = os.path.join(geodir,'coffs')
        logfile = os.path.join(procdir,'log',
                               'offset_fitm_{0}.log'.format(masterstr))
        print('Fitting offsets...')
        if not offset_fitm(offfile,ccpfile,diffpar,coffs,coffs+'ets','0.2','1',logfile):
            print("\nError:", file=sys.stderr)
            print("Something went wrong during the offset function fittin.", file=sys.stderr)
            return 3 
    ############################################################ Refine lookup table
        demwidth,demlength = get_dem_size(os.path.join(geodir,'EQA.dem_par')) # isn't this demseg?
        lutfine = os.path.join(geodir,masterstr+'.lt_fine')
        logfile = os.path.join(procdir,'log',
                               'gc_map_fine_{0}.log'.format(masterstr))
        if not gc_map_fine(lut,demwidth,diffpar,lutfine,'1',logfile):
            # Error in lookup table refining
            print("\nError:", file=sys.stderr)
            print("Something went wrong refining the lookup table.", file=sys.stderr)
            return 4
############################################################ create dem seg mli
                                                #by gecoding master mli to demseg
    else:
        print('skipping fitting of MLI to DEM to improve geocoding')
        lutfine = lut
        demwidth,demlength = get_dem_size(os.path.join(geodir,'EQA.dem_par'))
    
    width, length = get_mli_size(mli2+'.par')
    demsegmli = os.path.join(geodir,'EQA.{0}.slc.mli'.format(masterstr))
    logfile = os.path.join(procdir,'log',
                           'geocode_back_{0}.log'.format(masterstr))
    print('Geocoding...')
    if not geocode_back(mli2,width,lutfine,demsegmli,demwidth,
                                        demlength,'2','0',logfile):
        # Error in lookup table refining
        print("\nError:", file=sys.stderr)
        print("Something went wrong creating the DEM mli.", file=sys.stderr)
        return 5
############################################################ Create a height and 
                                    #intensisty ras of the mli2 data based on DEM
    logfile = os.path.join(procdir,'log',
                           'rashgt_{0}_hgt.log'.format(masterstr))
    if int(width)/1000 > 1:
        reducfac = str(int(width)/1000)
    else:
        reducfac = '1'
    if not rashgt(demseg,demsegmli,str(demwidth),'-','-','-',reducfac,reducfac,'500',logfile):
        print("\nError:", file=sys.stderr)
        print("Something went wrong creating the DEM mli sunraster file", file=sys.stderr)
        return 5
############################################################ Geocode (resample) DEM 
                                                            # -> mli2 coords
    hgtfile = os.path.join(geodir,masterstr+'.hgt')
    logfile = os.path.join(procdir,'log',
                           'geocode_{0}.log'.format(masterstr))
    if not geocode(lutfine,demseg,str(demwidth),hgtfile,str(width),str(length),'2','0',logfile):
        print("\nError:", file=sys.stderr)
        print("Something went wrong creating the master heightfile", file=sys.stderr)
        return 6
############################################################ Create a height and 
                                                            #intensisty ras of the
                                                            #mli2 data based on DEM
    if not rashgt(hgtfile,mli2,str(width),'-','-','-',reducfac,reducfac,'500',logfile):
        print("\nError:", file=sys.stderr)
        print("Something went wrong creating the DEM mli sunraster file", file=sys.stderr)
        return 6
    return 0


def get_nomissing_rslcs(rslcCache, mstrDate, builtRslcs):
    if type(builtRslcs) == type([]):
        builtRslcs = pd.DataFrame(builtRslcs)
    builtRslcs_nomissing = pd.DataFrame()
    try:
        masterstr = mstrDate.strftime('%Y%m%d')
    except:
        #in case it is already masterstr
        masterstr = mstrDate
    master_rslc = os.path.join(rslcCache, masterstr, masterstr+'.rslc')
    #if os.path.exists(master_rslc):
        #size_master = os.path.getsize(master_rslc)
    size_master_rslcs = 0
    for iwrslc in glob(os.path.join(rslcCache, masterstr, masterstr+'.IW?.rslc')):
        size_iwrslc = os.path.getsize(iwrslc)
        size_master_rslcs = size_master_rslcs + size_iwrslc
    for i,rslcdate in builtRslcs.iterrows():
        try:
            rslcdate_str = pd.Timestamp(rslcdate.values[0]).strftime('%Y%m%d')
        except:
            #to make it work for list of epoch strings
            rslcdate_str = rslcdate.values[0]
        rslcfile = os.path.join(rslcCache, rslcdate_str, rslcdate_str+'.rslc')
        #if os.path.exists(rslcfile):
        #    size_rslc = os.path.getsize(rslcfile)
        #    if size_rslc == size_master:
        #        builtRslcs_nomissing = builtRslcs_nomissing.append(rslcdate)
        size_rslcs = 0
        # let's check sizes if there is no mosaic - e.g. if we unzipped the file from LiCSAR_procdir
        for iwrslc in glob(os.path.join(rslcCache, rslcdate_str, rslcdate_str+'.IW?.rslc')):
            size_iwrslc = os.path.getsize(iwrslc)
            size_rslcs = size_rslcs + size_iwrslc
        if size_rslcs == size_master_rslcs:
            builtRslcs_nomissing = builtRslcs_nomissing.append(rslcdate)
    #else:
    #    print('ERROR - master RSLC mosaic does not exist!')
    return builtRslcs_nomissing


################################################################################
#Get DEM size function
################################################################################
def get_dem_size(demfile):
    """ Extracts the DEM width and length from parameter file
    """
    width = grep1('width:',demfile).split(':')[1].strip()
    length = grep1('nlines:',demfile).split(':')[1].strip()
    return [int(width),int(length)]

################################################################################
#Mosiac rslc
################################################################################
def mosaic_rslc(procdir,slavedate,masterdate,rglks,azlks,swathlist):
    """
    Wraps around SLC_mosaic_S1_TOPS - to mosaic subswathes into a single rslc
    """
    masterrslcdir = os.path.join(procdir,'RSLC',masterdate.strftime('%Y%m%d'))
    slaverslcdir = os.path.join(procdir,'RSLC',slavedate.strftime('%Y%m%d'))
    #master tab file
    masterslctab = os.path.join(procdir,'tab',masterdate.strftime('%Y%m%d')+'_tab')
    masterfilename = os.path.join(masterrslcdir,masterdate.strftime('%Y%m%d')+'.rslc')
    rc, msg = make_SLC_tab(masterslctab,masterfilename,swathlist)
    if rc > 0:
        print('Something went wrong creating the master tab file...')
        return 1
    #slave slc tab file
    slaverslctab = os.path.join(procdir,'tab',slavedate.strftime('%Y%m%d')+'_tab')
    slavefilename = os.path.join(slaverslcdir,slavedate.strftime('%Y%m%d')+'.rslc')
    rc, msg = make_SLC_tab(slaverslctab,slavefilename,swathlist)
    if rc > 0:
        print('Something went wrong creating the slave tab file...')
        return 1
    logfile = os.path.join(procdir,'log',"mosaic_rslc_{0}.log".format(slavedate.strftime('%Y%m%d')))
    if SLC_mosaic_S1_TOPS(slaverslctab,slavefilename,rglks,azlks,logfile,mastertab=masterslctab):
        # coreg was done ok, mosaicked ok, so geocode it now...
        rc = os.system('create_geoctiffs_to_pub.sh -M {0} {1}'.format(os.path.realpath(procdir), slavedate.strftime('%Y%m%d')))
        return 0
    else:
        return 1

################################################################################
#Rebuild mosaiced rslc
################################################################################
def rebuild_rslc(procdir,slavedate,masterdate,rglks,azlks):
    """
    Recreates (mosaics) a rslc from the subswathe rslc's iff the slave rslc doesn't
    exist. 
    """
    masterrslcdir = os.path.join(procdir,'RSLC',masterdate.strftime('%Y%m%d'))
    slaverslcdir = os.path.join(procdir,'RSLC',slavedate.strftime('%Y%m%d'))
    masterIWRSLC = glob(masterrslcdir+masterdate.strftime('/%Y%m%d.IW[1-3].rslc'))
    slaveIWRSLC = glob(slaverslcdir+slavedate.strftime('/%Y%m%d.IW[1-3].rslc'))
    slaverslc = os.path.join(slaverslcdir,slavedate.strftime('%Y%m%d.rslc'))
    if os.path.exists(slaverslc):
        print("slave rslc already exists")
        return 3
    if masterIWRSLC and slaveIWRSLC:
        swathlist = list()
        for siwr in slaveIWRSLC:
            mtch = re.search('.*(IW[1-3]).*',siwr)
            swathlist.append(mtch.groups()[0])
        swathlist.sort()
        if mosaic_rslc(procdir,slavedate,masterdate,rglks,azlks,swathlist):
            return 1
        else:
            print("Succesfully rebuilt rslc from subswath rslcs")
            return 0
    else:
        print("Could not find master or slave IW*.rslc, to remosaic")
        return 2

################################################################################
#get_mli_size function
################################################################################
def get_mli_size(mlifile):
    """ Extracts the DEM width and length from parameter file
    """
    width = grep1('range_samples:',mlifile).split(':')[1].strip()
    length = grep1('azimuth_lines:',mlifile).split(':')[1].strip()
    return [int(width),int(length)]

################################################################################
# Co-register slave function
################################################################################
#def coreg_slave_common(slavedate,slcdir,rslcdir,masterdate,framename,procdir, lq, job_id):
def coreg_slave_common(procdir,masterdate,masterrslcdir,slavedate,slaveslcdir,slaverslcdir,slaveslctab,slaverfilename,slave3tab,qualityfile,crop = False, eidp = False):
    # so far we were using same procedure as original GAMMA's S1_coreg_TOPS
    # the only introduced difference was using intensity cross correlation, but this seems not needed
    # so in order to keep the processing GAMMA-compatible also for future, let's use their original script
    # 
    # the original script is however updated - following changes were done:
    # - the generation of interferogram has been removed
    # - max 2 iterations for intensity matching (up to 0.01 px)
    # - SD correction is to 0.0005 px
    # but otherwise everything is kept same and saved to S1_coreg_TOPS_noifg
    # export OMP_NUM_THREADS=1
    if not crop:
        croptext=''
        geodir = os.path.join(procdir,'geo')
    else:
        croptext='_crop_'+slavedate.strftime('%Y%m%d')
        geodir = os.path.join(slaverslcdir,'geo')
    if slave3tab:
        slave3datestr = slave3tab.split('/')[-1][0:8]+croptext
    else:
        slave3datestr = ''
    #slaverslcdir = os.path.join(procdir,'RSLC',slavedate.strftime('%Y%m%d'))
    if not os.path.exists(slaverslcdir):
        os.mkdir(slaverslcdir)
    masterslctab = os.path.join(procdir,'tab',masterdate.strftime('%Y%m%d')+croptext+'_tab')
    #slaveslctab = os.path.join(procdir,'tab',slavedate.strftime('%Y%m%d')+'_tab')
    slaverslctab = os.path.join(procdir,'tab',slavedate.strftime('%Y%m%d')+'R_tab')
    demhgt = os.path.join(geodir,masterdate.strftime('%Y%m%d')+croptext+'.hgt')
    masterdatestr = masterdate.strftime('%Y%m%d')+croptext
    slavedatestr = slavedate.strftime('%Y%m%d')+croptext
    #pair = os.path.join(slaverslcdir,masterdate.strftime('%Y%m%d')+'_'+
    #                       slavedate.strftime('%Y%m%d'))
    logfile = os.path.join(procdir,'log','S1_coreg_TOPS_'+masterdate.strftime('%Y%m%d')+'_'+slavedate.strftime('%Y%m%d')+'.log')
    errfile = os.path.join(procdir,'log','S1_coreg_TOPS_'+masterdate.strftime('%Y%m%d')+'_'+slavedate.strftime('%Y%m%d')+'.err')
    
    
    #check if we have all IWs! edit tabs if not
    stab_pd = pd.read_csv(slaveslctab, header = None, sep=' ')
    
    iws = []
    for i, row in stab_pd.iterrows():
        if not os.path.exists(row[0]):
            print('skipping link to non-existing slc file: '+row[0])
            iws.append(os.path.basename(row[0]).split('.')[1])
    if iws:
        masterslctab2 = masterslctab+'.'+slavedate.strftime('%Y%m%d')+'.noiw'
        slaveslctab2 = slaveslctab+'.noiw'
        slaverslctab2 = slaverslctab+'.noiw'
        if slave3tab:
            slave3tab2 = slave3tab+'.noiw'
            shutil.copy(slave3tab,slave3tab2)
        shutil.copy(masterslctab,masterslctab2)
        shutil.copy(slaveslctab,slaveslctab2)
        shutil.copy(slaverslctab,slaverslctab2)
        for iw in iws:
            os.system("sed -i '/{}/d' {}".format(iw,masterslctab2))
            os.system("sed -i '/{}/d' {}".format(iw,slaveslctab2))
            os.system("sed -i '/{}/d' {}".format(iw,slaverslctab2))
            if slave3tab:
                os.system("sed -i '/{}/d' {}".format(iw,slave3tab2))
        masterslctab = masterslctab2
        slaveslctab = slaveslctab2
        slaverslctab = slaverslctab2
        if slave3tab:
            slave3tab = slave3tab2
        print('not using hgt here')
        demhgt = '-'
    # S1_coreg_TOPS <SLC1_tab> <SLC1_ID> <SLC2_tab> <SLC2_ID> <RSLC2_tab> [hgt] [RLK] [AZLK] [poly1] [poly2] [cc_thresh] [fraction_thresh] [ph_stdev_thresh] [cleaning] [flag1] [RSLC3_tab] [RSLC3_ID]
    #cmd = 'S1_coreg_TOPS_noifg {} {} {} {} {} {} {} {} - - 0.8 0.01 0.8 0 1 {} {} > {} 2> {}'.format(masterslctab, masterdatestr, slaveslctab, slavedatestr, 
    #            slaverslctab, demhgt, str(gc.rglks), str(gc.azlks), slave3tab, slave3datestr, logfile, errfile)
    #ML - gamma since version 20200728 - this really improved - we can coregister data we were not able to before!
    # usage: ScanSAR_coreg.py <SLC1_tab> <SLC1_ID> <SLC2_tab> <SLC2_ID> <RSLC2_tab> [hgt] [rlks] [azlks] [options]
    extracmd = ' '
    if slave3datestr:
        extracmd = '--RSLC3_tab {0} --RSLC3_ID {1}'.format(slave3tab, slave3datestr)
    #cmd = 'ScanSAR_coreg.py {0} {1} {2} {3} {4} {5} {6} {7} --use_existing --no_int --no_cleaning --wdir {8} '.format(masterslctab, masterdatestr, slaveslctab, slavedatestr, slaverslctab, demhgt, str(gc.rglks), str(gc.azlks), procdir) + extracmd + ' > {0} 2> {1}'.format(logfile, errfile)
    cmd = 'ScanSAR_coreg.py {0} {1} {2} {3} {4} {5} {6} {7} --no_int --wdir {8} '.format(masterslctab, masterdatestr, slaveslctab, slavedatestr, slaverslctab, demhgt, str(gc.rglks), str(gc.azlks), procdir) + extracmd + ' > {0} 2> {1}'.format(logfile, errfile)
    print(cmd)
    #now this may take some 80 minutes
    rc = os.system(cmd)
    # finally do the checks and move files where they belong
    gamma_qual = os.path.join(procdir, masterdatestr+'_'+slavedatestr+'.coreg_quality')
    gamma_lut = os.path.join(procdir, slavedatestr+'.mli.lt')
    gamma_off = os.path.join(procdir, masterdatestr+'_'+slavedatestr+'.off')
    gamma_coreg_results = os.path.join(procdir, masterdatestr+'_'+slavedatestr+'.results')
    pom = 1
    if os.path.exists(gamma_coreg_results):
        resfile = os.path.join(procdir,'RSLC',slavedate.strftime('%Y%m%d'),masterdate.strftime('%Y%m%d')+'_'+slavedate.strftime('%Y%m%d')+'.results')
        rc = shutil.copyfile(gamma_coreg_results, resfile)
    else:
        print('error - the results wile does not exist: '+gamma_coreg_results)
        #pom = 1
    #else:
    #    pom = 0
    if os.path.exists(gamma_off):
        offfile = os.path.join(procdir,'RSLC',slavedate.strftime('%Y%m%d'),masterdate.strftime('%Y%m%d')+'_'+slavedate.strftime('%Y%m%d')+'.off')
        rc = shutil.copyfile(gamma_off,offfile)
        #pom = 1
    else:
        pom = 0
    if os.path.exists(gamma_lut):
        lutfile = os.path.join(procdir,'RSLC',slavedate.strftime('%Y%m%d'),masterdate.strftime('%Y%m%d')+'_'+slavedate.strftime('%Y%m%d')+'.slc.mli.lt')
        rc = shutil.copyfile(gamma_lut,lutfile)
        #pom = 1
    else:
        print('output lut file was not generated. this happens often when SD is wrongly estimated')
        pom = 0
    if not os.path.exists(gamma_qual):
        pom = 0
    if pom == 0:
        print('Error - some of the result files were not generated')
        return 3
    # fill the coreg_quality file:
    with open(qualityfile, "a") as myfile:
            rc = myfile.write(str(cmd)+'\n')
            rc = myfile.write("Improvement of coregistration offset using intensity cross-correlation:\n")
    #cmd = 'grep matching_iteration_ {} >> {}'.format(gamma_qual, qualityfile)
    cmd = 'grep -A3 "matching iteration" {} >> {}'.format(gamma_qual, qualityfile)
    rc = os.system(cmd)
    with open(qualityfile, "a") as myfile:
            rc = myfile.write("Iterative improvement of refinement offset azimuth overlap regions:\n")
    #cmd = 'grep az_ovr_iteration_ {} >> {}'.format(gamma_qual, qualityfile)
    cmd = 'grep "az_ovr iteration" {} >> {}'.format(gamma_qual, qualityfile)
    rc = os.system(cmd)
    total_daztmp = grep1('azimuth_offset_polynomial:',gamma_off)
    total_daz = float(total_daztmp.split()[1])
    #if this is 0, it means an error in ESD estimation and should be cancelled - but keeping it 'bad' for NRT..
    if (total_daz == 0):
        print('daz was estimated as 0.0 px - this means for GAMMA that it failed J')
        if eidp:
            print('keeping it now, for the responder')
        else:
            return 2
    with open(qualityfile, "a") as myfile:
            rc = myfile.write("Total azimuth offset (w.r.t. LUT): {} (daz in SLC pixel)\n".format(total_daz))
            rc = myfile.write("(if the last iteration led to |daz| < 0.0005 px then this iteration was ignored)\n")
    return 0

#########################################################################

def coreg_slave(slavedate,slcdir,rslcdir,masterdate,framename,procdir, lq, job_id, maxdays_sd = 181, eidp = False):
    """ Coregister and resample slave to master geometry
    masterdate is of type dt.datetime.date...
    """
    print('\nCoregistering slave {0}...'.format(slavedate.date()))
    #get/create slave slc directory paths
    slaveslcdir = os.path.join(slcdir,slavedate.strftime('%Y%m%d'))
    slaverslcdir = os.path.join(rslcdir,slavedate.strftime('%Y%m%d'))
    if os.path.exists(slaverslcdir):
        shutil.rmtree(slaverslcdir)
    os.mkdir(slaverslcdir)
    #Create lock file for coregistration
    slaveLockFile = slaverslcdir+slavedate.strftime('/%Y%m%d.lock')
    open(slaveLockFile,'a').close() # OSX doesn't like os.mknod
    #master slc paths
    masterslcdir = os.path.join(slcdir,masterdate.strftime('%Y%m%d'))
    masterrslcdir = os.path.join(rslcdir,masterdate.strftime('%Y%m%d'))
    #calc diff date
    masterdiff = abs(slavedate.date()-masterdate)
    print("master difference is {0}".format(masterdiff))
    #create slave diff list
    procslavelist = [s for s in os.listdir(rslcdir) if 
            s != masterdate.strftime( '%Y%m%d') and s != slavedate.strftime('%Y%m%d')]
    procslavelist = get_nomissing_rslcs(rslcdir, masterdate, procslavelist)
    if not procslavelist.empty:
        procslavelist = procslavelist.values.transpose().tolist()[0]
    else:
        procslavelist = []
    print("Found existing RSLCs with probably no missing bursts: ... {0}".format(procslavelist))
    procslavediff = [dt.datetime.strptime(s,'%Y%m%d')-slavedate for s in
                                                                    procslavelist]
    # print "proc slave diffs {0}".format(procslavediff)
    #get slave bursts
    slavebursts = lq.get_frame_bursts_on_date(framename,slavedate)
    #get master bursts - should be same as for frame!
    masterbursts = lq.get_bursts_in_frame(framename)
    #try:
    #    masterbursts = lq.get_frame_bursts_on_date(framename,masterdate)
    #except:
    #    print('warning, cannot get list of master bursts. assuming same bursts as the secondary img..')
    #    masterbursts = slavebursts
    #    #i may do some more check here... like in arc_licsar.sh
############################################################ Find nearest (date) coregestered slave (already processed slave) to use as auxiliary
    rslc3override = ''
    cond = True
    while cond: #sweep through processed slaves till an appropriate aux is found, or till no processed slave available
        if procslavediff: # other slaves have been processed
            if np.min(np.abs(procslavediff)) < masterdiff and np.min(
                    np.abs(procslavediff)) > dt.timedelta(0): 
                # slaves closer than master
                slave3ix = np.argsort(np.abs(procslavediff))[0]
                slave3date = dt.datetime.strptime(procslavelist[slave3ix],'%Y%m%d')
                print("found potential aux slave date {0}".format(slave3date))
                if rslc3override == '':
                    rslc3override = slave3date
                slave3rslcdir = os.path.join(rslcdir,slave3date.strftime('%Y%m%d'))
                #print(slavebursts)
                #print(masterbursts)
                try:
                    slave3bursts = lq.get_frame_bursts_on_date(framename,slave3date)
                except:
                    #happens in EIDP when the -z filelist doesn't include the RSLC3
                    slave3bursts = ''
                if len(slave3bursts) == 0:
                    print('warning, cannot get list of RSLC3 bursts. assuming same bursts as master..')
                    print('THIS ACTUALLY IS CONSIDERED A BUG, SO CAREFUL!')
                    slave3bursts = masterbursts
                if os.path.exists(slave3rslcdir+slave3date.strftime('/%Y%m%d.lock')):
                        print("found lock {0}, not using date".format(slave3date))
                if glob(slave3rslcdir+'/*.IW[1-3].rslc') and not os.path.exists(slave3rslcdir+slave3date.strftime('/%Y%m%d.lock')):
                    missing = False
                    for sb in slavebursts:
                        if not sb in slave3bursts:
                            missing = True # aux slave does not have all of slave bursts
                else:
                    missing = True
                if missing: # Remove missing burst slave from list and try again
                    print("not valid as aux slave")
                    if len(procslavelist)>1:
                        procslavelist.pop(slave3ix)
                        procslavediff.pop(slave3ix)
                    else:
                        cond = False
                        print("but it is the last available - trying to keep this one, see what happens (inform Milan if it fails)")
                else: # Accept slave
                    print("using aux slave {0}".format(slave3date))
                    cond = False
            else: # No slaves closer than master
                slave3date = ''
                print("no aux slave used")
                cond = False
        else: # No slaves processed yet, or none left in list
            slave3date = '' #don't use aux slave.
            print("no aux slave used")
            cond = False
    #check to limit max btemp
    if not slave3date:
        tocheck = masterdate
    else:
        tocheck = slave3date.date()
    if abs(tocheck-slavedate.date()).days > maxdays_sd:
        print('time difference between RSLCs exceed max Btemp allowed for SD estimation = '+str(maxdays_sd)+' days')
        if not eidp:
            print('skipping')
            os.remove(slaveLockFile)
            return 1
        else:
            if not rslc3override == '':
                print('trying to override the bursts check - but expect problems..')
                slave3date = rslc3override
    
    #Get sorted list of swaths
    swathlist = [x[0].split('_')[1] for x in masterbursts]
    swathlist = set(swathlist)
    swathlist = list(swathlist)
    swathlist.sort()
    #Check for bursts missing in slaves but present in master
    missingbursts = []
    for mb in masterbursts:
        if not mb in slavebursts:
            missingbursts.append(mb)
############################################################ Create input tab files
    #master tab file
    print('WARNING, we will now not generate master tab file - it may contain wrong path..')
    masterslctab = os.path.join(procdir,'tab',masterdate.strftime('%Y%m%d')+'_tab')
    masterfilename = os.path.join(masterrslcdir,masterdate.strftime('%Y%m%d')+'.rslc')
    if not os.path.exists(masterslctab):
        rc, msg = make_SLC_tab(masterslctab,masterfilename,swathlist)
        if rc > 0:
            print('Something went wrong creating the master tab file...')
            return 1
    #slave slc tab file
    slaveslctab = os.path.join(procdir,'tab',slavedate.strftime('%Y%m%d')+'_tab')
    slavefilename = os.path.join(slaveslcdir,slavedate.strftime('%Y%m%d')+'.slc')
    if not os.path.exists(slaveslctab):
        rc, msg = make_SLC_tab(slaveslctab,slavefilename,swathlist)
        if rc > 0:
            print('Something went wrong creating the slave tab file...')
            return 1
    #resampled slave slc tab file
    slaverslctab = os.path.join(procdir,'tab',slavedate.strftime('%Y%m%d')+'R_tab')
    slaverfilename = os.path.join(slaverslcdir,slavedate.strftime('%Y%m%d')+'.rslc')
    if not os.path.exists(slaverslctab):
        rc, msg = make_SLC_tab(slaverslctab,slaverfilename,swathlist)
        if rc > 0:
            print('Something went wrong creating the slave resampled tab file...')
            return 1
    #auxillary slave tab file -> used for refinement passes
    if slave3date:
        slave3tab = os.path.join(procdir,'tab',slave3date.strftime('%Y%m%d_tab'))
        slave3filename = os.path.join(rslcdir,
                                          slave3date.strftime('%Y%m%d'),
                                          slave3date.strftime('%Y%m%d')+'.rslc')
        rc, msg = make_SLC_tab(slave3tab,slave3filename,swathlist)
        if rc > 0:
            print('Something went wrong creating the auxiliary slave tab file...')
            return 1
    else:
        slave3tab = ''
############################################################ Geometric coregistration
    qualityfile=os.path.join(procdir,'log','coreg_quality_'+
                              masterdate.strftime('%Y%m%d')+'_'+
                              slavedate.strftime('%Y%m%d')+'.log')
    with open(qualityfile, "w") as myfile:
        rc = myfile.write('Coregistration quality of SLC dated '+slavedate.strftime('%Y%m%d')+"\n")
    if slave3date:
        with open(qualityfile, "a") as myfile:
            rc = myfile.write("Spectral diversity estimation will be performed through near RSLC: "+slave3date.strftime('%Y%m%d')+"\n")
    #print('WARNING - some issues detected with wrong detection of missing bursts')
    #print('doing workaround here - ALWAYS check the burst correspondence through GAMMA command')
    #print('this needs to be improved - I recommend full rewrite of coreg_lib.py -- really, better do it from scratch as this contains routines from year 0..er.. 2015')
    #if not missingbursts:
    #    missingbursts = ['fakemissingburst']
    #if not missingbursts:
    if 1==1:        #print('All bursts available, no recropping of master necessary...')
        print('debug place 1')
        #coreg_slave_common(slavedate,slcdir,rslcdir,masterdate,framename,procdir, lq, job_id)
        rc = coreg_slave_common(procdir,masterdate,masterrslcdir,slavedate,slaveslcdir,slaverslcdir,slaveslctab,slaverfilename,slave3tab,qualityfile, eidp = eidp)
        if rc != 0:
            print("\nError:", file=sys.stderr)
            print("Something went wrong during coregistration", file=sys.stderr)
            if os.path.exists(slaverslcdir):
                shutil.rmtree(slaverslcdir)
            return rc
        #Do the RSLC mosaic (should actually exist from the previous steps..
        if not os.path.exists(slaverfilename):
            logfilename = os.path.join(procdir,'log','SLC_mosaic_S1_TOPS_'+
                                   slavedate.strftime('%Y%m%d')+
                                   '.log')
            if not SLC_mosaic_S1_TOPS(slaverslctab,slaverfilename,gc.rglks,gc.azlks,logfilename,mastertab=masterslctab):
                print('Something went wrong mosaicing subswaths '\
                    'together. Log file {0}. Continuing with next acquisition '\
                    'date.'.format(logfilename), file=sys.stderr)
                return 3
        #Multilook (average) resampled slc
        logfilename =  os.path.join(procdir,'log',
                                            'multilookRSLC_{0}_crop.log'.format(
                                                slavedate.strftime('%Y%m%d')))
        multicall = 'multilookRSLC {0} {1} {2} 1 {3} &> {4}'.format(
                slavedate.strftime('%Y%m%d'),gc.rglks,gc.azlks,slaverslcdir,logfilename)
        rc = os.system(multicall)
        if rc != 0:
            print('Something went wrong multilooking the '\
                    'resampled slave image. Log file {0}. Continuing with '\
                    'next acquisition date.'.format(logfilename), file=sys.stderr)
            return rc
                
    else:
############################################################ Crop master to fit with smaller slave
        # this solution was working also with the 20181130 gamma codes. Keeping as it is..
        #print('There are '+str(len(missingbursts))+' missing bursts at either or both ends of the scene')
        #print('..recropping master file')
        #if slave3date:
        #    print('..and also recropping auxilliary slave')
        #with open(qualityfile, "a") as myfile:
        #    rc = myfile.write("There were missing bursts. Related coreg routine was used \n")
        masterslctab = os.path.join(procdir,'tab',masterdate.strftime('%Y%m%d')
                                                                            +'_tab')
        masterslcfilename = os.path.join(masterrslcdir,masterdate.strftime('%Y%m%d')
                                                                            +'.rslc')
        croptext = '_crop_'+slavedate.strftime('%Y%m%d')
        croptab = os.path.join(procdir,'tab',masterdate.strftime('%Y%m%d')+croptext+
                                                                        '_tab')
        cropfilename = os.path.join(masterrslcdir,masterdate.strftime('%Y%m%d')
                                                                        +croptext+'.rslc')
        croptext = '_crop_'+slavedate.strftime('%Y%m%d')
        slave_croptab = os.path.join(procdir,'tab',slavedate.strftime('%Y%m%d')+croptext+
                                                                        '_tab')
        slave_cropfilename = os.path.join(slaverslcdir,
                                slavedate.strftime('%Y%m%d')+croptext+'.rslc')
        # prepare tab:
        rc, msg = make_SLC_tab(slave_croptab,slave_cropfilename,swathlist)
        if not os.path.exists(slaverslcdir):
            os.mkdir(slaverslcdir)
        #print('..rebuilding slave file')
        # if the master TOPS_par is too old, need to regenerate it!
        #if .........:
        #     par_S1_SLC........
        # S1_coreg_TOPS_burst_selection <SLC1_tab> <SLC2_tab> <SLC2_tab1> [mode]
        #try new version:
        cmd = 'S1_coreg_TOPS_burst_selection {} {} {} 0 >/dev/null 2>/dev/null'.format(masterslctab, slaveslctab, slave_croptab)
        print('bursts selection based on GAMMA (may not work for old primary SLCs)')
        print(cmd)
        rc = os.system(cmd)
        if os.listdir(slaverslcdir):
            #so this went ok... now let's just use it for processing
            rc = coreg_slave_common(procdir,masterdate,masterrslcdir,slavedate,slaverslcdir,slaverslcdir,slave_croptab,slaverfilename,slave3tab,qualityfile, eidp=eidp)
            if rc != 0:
                print("\nError:", file=sys.stderr)
                print("Something went wrong during coregistration", file=sys.stderr)
                if os.path.exists(slaverslcdir):
                    shutil.rmtree(slaverslcdir)
                return rc
            #Do the RSLC mosaic (should actually exist from the previous steps..
            if not os.path.exists(slaverfilename):
                logfilename = os.path.join(procdir,'log','SLC_mosaic_S1_TOPS_'+
                                       slavedate.strftime('%Y%m%d')+
                                       '.log')
                if not SLC_mosaic_S1_TOPS(slaverslctab,slaverfilename,gc.rglks,gc.azlks,logfilename,mastertab=masterslctab):
                    print('Something went wrong mosaicing subswaths '\
                        'together. Log file {0}. Continuing with next acquisition '\
                        'date.'.format(logfilename), file=sys.stderr)
                    return 3
            #Multilook (average) resampled slc
            logfilename =  os.path.join(procdir,'log',
                                                'multilookRSLC_{0}_crop.log'.format(
                                                    slavedate.strftime('%Y%m%d')))
            multicall = 'multilookRSLC {0} {1} {2} 1 {3} &> {4}'.format(
                    slavedate.strftime('%Y%m%d'),gc.rglks,gc.azlks,slaverslcdir,logfilename)
            rc = os.system(multicall)
            if rc != 0:
                print('Something went wrong multilooking the '\
                        'resampled slave image. Log file {0}. Continuing with '\
                        'next acquisition date.'.format(logfilename), file=sys.stderr)
                return rc
        else:
            #this often means that the master file was initialised using old version of GAMMA
            #loop through swaths
            print('(cropping on swath level)')
            for iw in swathlist: #['IW1','IW2','IW3']:
            #get master and slave bursts
                iwthisburst_m = sorted([b[0] for b in masterbursts if iw in b[0]])
                iwthisburst_s = sorted([b[0] for b in slavebursts if iw in b[0]])
            #find indices which match slave bursts
                iwstart = iwthisburst_m.index(iwthisburst_s[0])+1
                iwstop = iwthisburst_m.index(iwthisburst_s[-1])+1
            #Create master slc tab file
                mastertabiw = masterslctab+'.'+iw
                rc, msg = make_SLC_tab(mastertabiw,masterslcfilename,[iw])
                if rc > 0:
                    print('Something went wrong creating the master tab file...')
                    return 1
            #create a corresponding cropped master slc tab file
                croptabiw = croptab+'.'+iw
                rc, msg = make_SLC_tab(croptabiw,cropfilename,[iw])
                if rc > 0:
                    print('Something went wrong creating the cropped master tab file...')
                    return 1
                bursttab = os.path.join( procdir,'tab', 'master'+croptext+'_{0}_burst_tab'.format( iw ) )
                rc, msg = make_burst_tab( bursttab, iwstart, iwstop)
                if rc > 0:
                    print('\nProblem creating burst '\
                                        'tab {0} for date {1}. Error '\
                                        'message:'\
                                        '\n {2}'\
                                        '\nContinuing with next '\
                                        'date.'.format( bursttab, 
                                                date, msg ), file=sys.stderr)
                    return 2
                logfilename = os.path.join(procdir,'log','SLC_copy_S1_TOPS_'+
                                           masterdate.strftime('%Y%m%d')+
                                           croptext+'.log')
            #copy out bursts which are relevant to the cropped master
                if SLC_copy_S1_TOPS(mastertabiw,croptabiw, bursttab,'',procdir,
                                                                        logfilename):
                    print("\nError:", file=sys.stderr)
                    print("Something went wrong recropping the master SLC", file=sys.stderr)
                    return 5
            #if auxiliary slave is present, crop master to it as well
                slave3_croptab = ''
                if slave3date:
                    iwthisburst_s3 = sorted([b[0] for b in slave3bursts if iw in b[0]])
                    iwthisburst_s = sorted([b[0] for b in slavebursts if iw in b[0]])
            #find burst indices which correspond to slave limits
                    iwstart = iwthisburst_s3.index(iwthisburst_s[0])+1
                    iwstop = iwthisburst_s3.index(iwthisburst_s[-1])+1
            #create aux slave tab file
                    slave3slctab = os.path.join(procdir,'tab',slave3date.strftime('%Y%m%d')
                                                                                    +'_tab')
                    slave3slcfilename = os.path.join(slave3rslcdir,
                                        slave3date.strftime('%Y%m%d')+'.rslc')
                    slave3slctabiw = slave3slctab+'.'+iw
                    rc, msg = make_SLC_tab(slave3slctabiw,slave3slcfilename,[iw])
                    if rc > 0:
                        print('Something went wrong creating the auxiliary slave tab file...')
                        return 1
                #create a corresponding crop tab file for aux slave
                    slave3_croptab = os.path.join(procdir,'tab',
                                    slave3date.strftime('%Y%m%d')+croptext+'_tab')
                    slave3_cropfilename = os.path.join(slave3rslcdir,
                                    slave3date.strftime('%Y%m%d')+croptext+'.rslc')
                    slave3_croptabiw = slave3_croptab+'.'+iw
                    rc, msg = make_SLC_tab(slave3_croptabiw,slave3_cropfilename,[iw])
                    if rc > 0:
                        print('Something went wrong creating the cropped '\
                                            'auxiliary slave tab file...')
                        return 1
                    logfilename = os.path.join(procdir,'log','SLC_copy_S1_TOPS_'+
                                               slave3date.strftime('%Y%m%d')+
                                               croptext+'.log')
                    bursttab = os.path.join( procdir,'tab', 'master_aux'+croptext+'_{0}_burst_tab'.format( iw ) )
                    rc, msg = make_burst_tab( bursttab, iwstart, iwstop)
                    if rc > 0:
                        print('\nProblem creating burst '\
                                            'tab {0} for date {1}. Error '\
                                            'message:'\
                                            '\n {2}'\
                                            '\nContinuing with next '\
                                            'date.'.format( bursttab, 
                                                    date, msg ), file=sys.stderr)
                        return 2
            #copy out relevent bursts in master to aux slavetab??? is it correct (Nick)??
            #       if SLC_copy_S1_TOPS(masterslctab,slave3_croptaw,.......
            #oh no, it was an error -- took too long to find it out..
                    if SLC_copy_S1_TOPS(slave3slctabiw,slave3_croptabiw,
                                       bursttab,'',procdir,logfilename):
                        print("\nError:", file=sys.stderr)
                        print("Something went wrong recropping the master SLC", file=sys.stderr)
                        return 5
            #mosaic cropped files
            logfilename = os.path.join(procdir,'log','SLC_mosaic_S1_TOPS_'+
                                       masterdate.strftime('%Y%m%d')+
                                       '_{0}'+croptext+'.log'.format(iw))
            rc, msg = make_SLC_tab(croptab,cropfilename,swathlist)
            print('(mosaicking cropped swaths)')
    
            if not SLC_mosaic_S1_TOPS(croptab,cropfilename,gc.rglks,gc.azlks,logfilename):
                print('Something went wrong mosaicing subswaths '\
                        'together. Log file {0}. Continuing with next acquisition '\
                        'date.'.format(logfilename), file=sys.stderr)
                return 3
            if slave3date:
                rc, msg = make_SLC_tab(slave3_croptab,slave3_cropfilename,swathlist)
                if not SLC_mosaic_S1_TOPS(slave3_croptab,slave3_cropfilename,gc.rglks,gc.azlks,logfilename):
                    print('Something went wrong mosaicing subswaths '\
                            'together. Log file {0}. Continuing with next acquisition '\
                            'date.'.format(logfilename), file=sys.stderr)
                    return 3
            #multilook (average) cropped data
            logfilename =  os.path.join(procdir,'log',
                                'multilookRSLC_'+masterdate.strftime('%Y%m%d')+croptext+'.log')
            print('(multilooking cropped mosaic)')
            multicall = 'multilookRSLC {0} {1} {2} 1 {3} &> {4}'.format(
                    masterdate.strftime('%Y%m%d')+croptext,gc.rglks,gc.azlks,
                    masterrslcdir,logfilename)
            rc = os.system(multicall)
            if rc != 0:
                print('Something went wrong multilooking the '\
                        'merged image. Log file {0}. Continuing with next '\
                        'acquisition date.'.format(logfilename), file=sys.stderr)
        #create file paths
            geocropdir = os.path.join(slaverslcdir,'geo')
            os.mkdir(geocropdir)
            demdir = os.path.join(procdir,'DEM')
            print('(cropping DEM, geocoding,..)')
            rc = geocode_dem(masterrslcdir,geocropdir,demdir,
                    procdir,masterdate.strftime('%Y%m%d')+croptext,gc.outres)
        #finally the coregistration itself
            rc = coreg_slave_common(procdir,masterdate,masterrslcdir,slavedate,slaveslcdir,slaverslcdir,slaveslctab,slaverfilename,slave3_croptab,qualityfile,True, eidp = eidp)
            if rc != 0:
                print("\nError:", file=sys.stderr)
                print("Something went wrong during coregistration", file=sys.stderr)
                if os.path.exists(slaverslcdir):
                    shutil.rmtree(slaverslcdir)
                return rc
        # as final touches, rslcs for all swaths separately will be padded. therefore there will be no issues mosaicking in the future
            for iw in swathlist: #['IW1','IW2','IW3']:
                if os.path.exists(os.path.join(slaverslcdir,slavedate.strftime('%Y%m%d')+'.'+iw+'.rslc')):
                    masterOrigRslc = os.path.join(masterrslcdir,masterdate.strftime('%Y%m%d')+'.'+iw+'.rslc')
                    masterCropRslc = os.path.join(masterrslcdir,masterdate.strftime('%Y%m%d')+croptext+'.'+iw +'.rslc')
                    if os.path.getsize(masterOrigRslc) != os.path.getsize(masterCropRslc):
                        # do the padding here... this seems like the easier solution on the 'ghosting' problem
                        # better (disk saving) option would be to save just azimuth offset (lines) towards orig. master
                        # and use it for generating RSLC mosaic. however this would need more coding and...
                        # we are not going to store the full RSLCs anyway....
                        # also, padding the IW RSLCs would be more savvy towards RAM memory needs
                        azoffsetcall = ['master2mastercrop_offset.sh',
                                        os.path.join(masterrslcdir,
                                                     masterdate.strftime('%Y%m%d')+'.'+iw),
                                        'rslc',
                                        os.path.join(masterrslcdir,
                                                     masterdate.strftime('%Y%m%d'))+croptext+'.'+iw,
                                        'rslc']
                        azoffset = subp.check_output(azoffsetcall).strip()
                        azoffset = azoffset.decode('ascii')
                        rslccall = ['rslc2rslc.sh',
                                    os.path.join(masterrslcdir,
                                                     masterdate.strftime('%Y%m%d')+'.'+iw),
                                    os.path.join(slaverslcdir,
                                                     slavedate.strftime('%Y%m%d')+'.'+iw),
                                    str(azoffset)]
                        rc = subp.check_call(rslccall)
                        if rc > 0:
                            print("\nError:", file=sys.stderr)
                            print("Something went wrong zero padding cropped RSLC.", file=sys.stderr)
                            return 4
            #Do the RSLC mosaic
            rc, msg = make_SLC_tab(masterslctab,masterslcfilename,swathlist)
            logfilename = os.path.join(procdir,'log','SLC_mosaic_S1_TOPS_'+
                                       slavedate.strftime('%Y%m%d')+
                                       '.log')
            if not SLC_mosaic_S1_TOPS(slaverslctab,slaverfilename,gc.rglks,gc.azlks,logfilename,mastertab=masterslctab):
                print('Something went wrong mosaicing subswaths '\
                        'together. Log file {0}. Continuing with next acquisition '\
                        'date.'.format(logfilename), file=sys.stderr)
                return 3
            #Multilook (average) resampled slc
            logfilename =  os.path.join(procdir,'log',
                                                'multilookRSLC_{0}_crop.log'.format(
                                                    slavedate.strftime('%Y%m%d')))
            multicall = 'multilookRSLC {0} {1} {2} 1 {3} &> {4}'.format(
                    slavedate.strftime('%Y%m%d'),gc.rglks,gc.azlks,slaverslcdir,logfilename)
            rc = os.system(multicall)
            if rc != 0:
                print('Something went wrong multilooking the '\
                        'resampled slave image. Log file {0}. Continuing with '\
                        'next acquisition date.'.format(logfilename), file=sys.stderr)
                return rc
    ####
    ############################################################ Pad cropped data
    #### this part is to save the azimuth offset (lines no.) between orig and cropped master. it may be useful for the LUT-only recoregistration
            #print('Padding coregistered data to the original master extent')
            azoffsetcall = ['master2mastercrop_offset.sh',
                            os.path.join(masterrslcdir,
                                         masterdate.strftime('%Y%m%d')),
                            'rslc',
                            os.path.join(masterrslcdir,
                                         masterdate.strftime('%Y%m%d'))+croptext,
                            'rslc']
            azoffset = subp.check_output(azoffsetcall).strip()
            azoffset = azoffset.decode('ascii')
            with open(qualityfile, "a") as myfile:
                myfile.write("There were missing bursts. The azimuth offset towards master SLC is: \n")
                myfile.write("azoffset_lines: "+str(azoffset)+"\n")
# this is residual from previous version where missing bursts were mosaicked and only then padded (here we would need to keep mosaic instead of IWs)
        #rslccall = ['rslc2rslc.sh',
        #            os.path.join(masterrslcdir,
        #                             masterdate.strftime('%Y%m%d')),
        #            os.path.join(slaverslcdir,
        #                             slavedate.strftime('%Y%m%d')),
        #            str(azoffset)]
        #rc = subp.check_call(rslccall)
        #if rc > 0:
        #    print("\nError:", file=sys.stderr)
        #    print("Something went wrong zero padding cropped RSLC.", file=sys.stderr)
        #    return 4
############################################################ remove lock file
    os.remove(slaveLockFile)

############################################################ Update job database (if not in batch mode)
    # Populate the RSLC table with the location of the new slave RSLC
    print("RSLC succesfully made, now to write to RSLC table...")
    print("job_id: %d" % job_id)


    if job_id != -1:
        print('Attempting to input new rslc SLAVE product')
        lq.set_new_rslc_product(job_id,slavedate.strftime('%Y%m%d'),
                                masterrslcdir+'/'+masterdate.strftime('%Y%m%d')+'.rslc',
                                slavedate.strftime('%Y%m%d') + '.rslc', slaverslcdir)

    return 0

################################################################################
# Re-Co-register slave function
################################################################################
def recoreg_slave(slavedate,slcdir,rslcdir,masterdate,framename,procdir,lq):
    """ Recoregister and resample slave to master geometry. This assumes lookup
    tables (and off files) are stored. keeping same parameters as the old version, for compatibility
    We expect lookup table to be in LUT/SLAVEDATE/MASTER_SLAVEDATE.slc.mli.lt (and off file should exist)
    Note the lt is after both rdc_trans+ICC coregistration, and the ESD estimate is inside the off file (in the same LUT/DATE directory)
    """
    print('\nRecoregistering epoch {0}...'.format(slavedate.date()))
    
    #get/create slave slc directory paths
    slaveslcdir = os.path.join(slcdir,slavedate.strftime('%Y%m%d'))
    slaverslcdir = os.path.join(rslcdir,slavedate.strftime('%Y%m%d'))
    slavelutdir = os.path.join(rslcdir,'../LUT',slavedate.strftime('%Y%m%d'))
    if not os.path.exists(slaverslcdir):
        os.mkdir(slaverslcdir)
        #return 6 # slave rslc dir doesn't exist!!

    #Create lock file for coregistration
    slaveLockFile = slaverslcdir+slavedate.strftime('/%Y%m%d.lock')
    open(slaveLockFile,'a').close() # OSX doesn't like os.mknod

    #master slc paths
    masterslcdir = os.path.join(slcdir,masterdate.strftime('%Y%m%d'))
    masterrslcdir = os.path.join(rslcdir,masterdate.strftime('%Y%m%d'))

    #this part may be needed to update! (?)
    #get slave bursts
    slavebursts = lq.get_frame_bursts_on_date(framename,slavedate)
    #get master bursts
    masterbursts = lq.get_frame_bursts_on_date(framename,masterdate)
    if not masterbursts:
        masterbursts = lq.get_bursts_in_frame(framename)
    #Get sorted list of swaths
    swathlist = [x[0].split('_')[1] for x in masterbursts]
    swathlist = set(swathlist)
    swathlist = list(swathlist)
    swathlist.sort()

    #Check for bursts missing in slaves but present in master
    missingbursts = []
    for mb in masterbursts:
        if not mb in slavebursts:
            missingbursts.append(mb)


############################################################ Create input tab files
    #master tab file
    masterslctab = os.path.join(procdir,'tab',masterdate.strftime('%Y%m%d')+'_tab')
    masterfilename = os.path.join(masterrslcdir,masterdate.strftime('%Y%m%d')+'.rslc')
    rc, msg = make_SLC_tab(masterslctab,masterfilename,swathlist)
    if rc > 0:
        print('Something went wrong creating the master tab file...')
        return 1


    #slave slc tab file
    slaveslctab = os.path.join(procdir,'tab',slavedate.strftime('%Y%m%d')+'_tab')
    slavefilename = os.path.join(slaveslcdir,slavedate.strftime('%Y%m%d')+'.slc')
    rc, msg = make_SLC_tab(slaveslctab,slavefilename,swathlist)
    if rc > 0:
        print('Something went wrong creating the slave tab file...')
        return 1

    #resampled slave slc tab file
    slaverslctab = os.path.join(procdir,'tab',slavedate.strftime('%Y%m%d')+'R_tab')
    slaverfilename = os.path.join(slaverslcdir,slavedate.strftime('%Y%m%d')+'.rslc')
    rc, msg = make_SLC_tab(slaverslctab,slaverfilename,swathlist)
    if rc > 0:
        print('Something went wrong creating the slave resampled tab file...')
        return 1

    #pair = os.path.join(slaverslcdir,masterdate.strftime('%Y%m%d')+'_'+
    #                       slavedate.strftime('%Y%m%d'))
    #offset parameter file path - should be the final one
    #offfile = pair+'.off'
    #refine2file = offfile+'.refine2'
    #refine2file = offfile
############################################################ Geomatric coregistration
    #if no missing bursts....
    #missingbursts = False
    if not missingbursts:
        #print('debug place 2')
        print('All bursts available, no recropping of master necessary...')
    #master mli param file path
        mastermlipar = os.path.join(masterslcdir,
                                 masterdate.strftime('%Y%m%d')+'.slc.mli.par')
    #master param file path
        masterpar = os.path.join(masterslcdir,
                                 masterdate.strftime('%Y%m%d')+'.slc.par')
    #slave mli param file path
        slavemlipar = os.path.join(slaveslcdir,
                                slavedate.strftime('%Y%m%d')+'.slc.mli.par')
    #master param file path
        slavepar = os.path.join(slaveslcdir,
                                slavedate.strftime('%Y%m%d')+'.slc.par')
    #DEM height file
        demhgt = os.path.join(procdir,'geo',masterdate.strftime('%Y%m%d')+'.hgt')
    #lookup table path
        lut = os.path.join(slavelutdir,
                           masterdate.strftime('%Y%m%d')+'_'+
                           slavedate.strftime('%Y%m%d')+'.slc.mli.lt')
        offfile = os.path.join(slavelutdir, masterdate.strftime('%Y%m%d')+'_'+
                           slavedate.strftime('%Y%m%d')+'.off')
        if not os.path.exists(lut):
            lut =  os.path.join(slaverslcdir,
                           masterdate.strftime('%Y%m%d')+'_'+
                           slavedate.strftime('%Y%m%d')+'.slc.mli.lt')
            offfile = os.path.join(slaverslcdir, masterdate.strftime('%Y%m%d')+'_'+
                           slavedate.strftime('%Y%m%d')+'.off')
            if not os.path.exists(lut):
                return 7 # Couldn't find previous look up table
        
        if not os.path.exists(offfile):
            licsar_procdir = os.environ['LiCSAR_procdir']
            track = str(int(framename[:3]))
            coreg_qual_old = os.path.join(licsar_procdir,track,framename,'log/coreg_quality_'+masterdate.strftime('%Y%m%d')+'_'+
                           slavedate.strftime('%Y%m%d')+'.log')
            if os.path.exists(coreg_qual_old):
                rc = regenerate_offfile(offfile, coreg_qual_old, mastermlipar)
            else:
                print('RSLC cannot be regenerated since neither off nor coreg_qual log file exists')
                return 7
            if not os.path.exists(offfile):
                return 7
        logfile = os.path.join(procdir,'log','rdc_trans_'+
                               masterdate.strftime('%Y%m%d')+'_'+
                               slavedate.strftime('%Y%m%d')+'.log')

        print('Resampling the SLC using the existing LUT')
        if not SLC_interp_lt_S1_TOPS(slaveslctab,slavepar,masterslctab,
                                     masterpar,lut,mastermlipar,slavemlipar,
                                     offfile,slaverslctab,slaverfilename,
                                     slaverfilename+'.par',logfile):
            print("\nError:", file=sys.stderr)
            print("Something went wrong resampling the slave SLC", file=sys.stderr)
            return 3
        #Multilook (average) resampled slc
        logfilename =  os.path.join(procdir,'log',
                                            'multilookRSLC_{0}_crop.log'.format(
                                                slavedate.strftime('%Y%m%d')))
        multicall = 'multilookRSLC {0} {1} {2} 1 {3} &> {4}'.format(
                slavedate.strftime('%Y%m%d'),gc.rglks,gc.azlks,slaverslcdir,logfilename)
        rc = os.system(multicall)
        if rc != 0:
            print('Something went wrong multilooking the '\
                    'resampled slave image. Log file {0}. Continuing with '\
                    'next acquisition date.'.format(logfilename), file=sys.stderr)

    else:
        mastermlipar = os.path.join(masterrslcdir,
                                 masterdate.strftime('%Y%m%d')+'_crop.rslc.mli.par')
        masterpar = os.path.join(masterrslcdir,
                                 masterdate.strftime('%Y%m%d')+'_crop.slc.par')
        slavemlipar = os.path.join(slaveslcdir,
                                slavedate.strftime('%Y%m%d')+'.slc.mli.par')
        slavepar = os.path.join(slaveslcdir,
                                slavedate.strftime('%Y%m%d')+'.slc.par')
        demhgt = os.path.join(slaverslcdir,'geo',masterdate.strftime('%Y%m%d')
                                                                +'_crop.hgt')
        lut = os.path.join(slaverslcdir,
                           masterdate.strftime('%Y%m%d')+'_'+
                           slavedate.strftime('%Y%m%d')+'.slc.mli.lt')
        offfile = os.path.join(slaverslcdir, masterdate.strftime('%Y%m%d')+'_'+
                           slavedate.strftime('%Y%m%d')+'.off')
        if not os.path.exists(lut):
            lut = os.path.join(slavelutdir,
                           masterdate.strftime('%Y%m%d')+'_'+
                           slavedate.strftime('%Y%m%d')+'.slc.mli.lt')
            offfile = os.path.join(slavelutdir, masterdate.strftime('%Y%m%d')+'_'+
                           slavedate.strftime('%Y%m%d')+'.off')
            if not os.path.exists(lut):
                return 7 # Couldn't find previous look up table
        ww=int(grep1('interferogram_width', offfile).split(':')[1].strip())
        ll=int(grep1('interferogram_azimuth_lines', offfile).split(':')[1].strip())
        dimensions_lut = ww*ll
        ww=int(grep1('range_samples', slavemlipar).split(':')[1].strip())
        ll=int(grep1('azimuth_lines', slavemlipar).split(':')[1].strip())
        dimensions_slavemli = ww*ll
        if not (dimensions_lut == dimensions_slavemli):
            print('This SLC has missing bursts not corresponding to the LUT (it was perhaps created from the full). Trying the full coreg here but better re-download the epoch data..')
            return coreg_slave(slavedate,slcdir,rslcdir,masterdate,framename,procdir, lq, job_id = -1, maxdays_sd = 181, eidp = False)
############################################################ Crop master to fit with smaller slave
        print('Missing bursts at either or both ends of the scene, recropping '\
                'master, and auxiliary slave if necessary...')
    #loop through swaths
        for iw in ['IW1','IW2','IW3']:
        #get master and slave bursts
            iwthisburst_m = sorted([b[0] for b in masterbursts if iw in b[0]])
            iwthisburst_s = sorted([b[0] for b in slavebursts if iw in b[0]])
        #find indices which match slave bursts
            i=0
            iwstart=False
            iwstop=False
            # Note: the '+1' here is because then we will extract bursts from GAMMA par files that number bursts starting from 1
            while i<len(iwthisburst_s):
                try:
                    iwstart = iwthisburst_m.index(iwthisburst_s[i])+1 
                    break
                except:
                    i=i+1
            i=0
            while i<len(iwthisburst_s):
                try:
                    iwstop = iwthisburst_m.index(iwthisburst_s[int(-1*(i+1))])+1
                    break
                except:
                    i=i+1
            if (not iwstart) or (not iwstop):
                print('common bursts not found, cannot recoregister - jumping to common coreg (or should we return error here?)')
                return coreg_slave(slavedate,slcdir,rslcdir,masterdate,framename,procdir, lq, job_id = -1, maxdays_sd = 181, eidp = False)
        #Create master slc tab file
            masterslctab = os.path.join(procdir,'tab',masterdate.strftime('%Y%m%d')
                                                                            +'_tab')
            masterslcfilename = os.path.join(masterrslcdir,masterdate.strftime('%Y%m%d')
                                                                            +'.rslc')
            rc, msg = make_SLC_tab(masterslctab,masterslcfilename,[iw])
            if rc > 0:
                print('Something went wrong creating the master tab file...')
                return 1
        #create a corresponding cropped master slc tab file
            croptab = os.path.join(procdir,'tab',masterdate.strftime('%Y%m%d')+
                                                                        'crop_tab')
            cropfilename = os.path.join(masterrslcdir,masterdate.strftime('%Y%m%d')
                                                                        +'_crop.rslc')
            rc, msg = make_SLC_tab(croptab,cropfilename,[iw])
            if rc > 0:
                print('Something went wrong creating the cropped master tab file...')
                return 1

            bursttab = os.path.join( procdir,'tab', 'master_crop_{0}_burst_tab'.format( iw ) )
            rc, msg = make_burst_tab( bursttab, iwstart, iwstop)
            if rc > 0:
                print('\nProblem creating burst '\
                                    'tab {0} for date {1}. Error '\
                                    'message:'\
                                    '\n {2}'\
                                    '\nContinuing with next '\
                                    'date.'.format( bursttab, 
                                            date, msg ), file=sys.stderr)
                return 2
            logfilename = os.path.join(procdir,'log','SLC_copy_S1_TOPS_'+
                                       masterdate.strftime('%Y%m%d')+
                                       '_crop.log')
        #copy out bursts which are relevant to the cropped master
            if SLC_copy_S1_TOPS(masterslctab,croptab,bursttab,'',procdir,
                                                                    logfilename):
                print("\nError:", file=sys.stderr)
                print("Something went wrong recropping the master SLC", file=sys.stderr)
                return 5

        logfilename = os.path.join(procdir,'log','SLC_mosaic_S1_TOPS_'+
                                   masterdate.strftime('%Y%m%d')+
                                   '_{0}_crop.log'.format(iw))

    #mosaic cropped file (does this mean if aux slave is present we only use that cropped file?)
        rc, msg = make_SLC_tab(croptab,cropfilename,swathlist)
        if not SLC_mosaic_S1_TOPS(croptab,cropfilename,gc.rglks,gc.azlks,logfilename):
            print('Something went wrong mosaicing subswaths '\
                    'together. Log file {0}. Continuing with next acquisition '\
                    'date.'.format(logfilename), file=sys.stderr)
            return 3
        logfilename =  os.path.join(procdir,'log',
                            'multilookRSLC_{0}_crop.log'.format(
                                masterdate.strftime('%Y%m%d')))
    #multilook (average) cropped data
        multicall = 'multilookRSLC {0} {1} {2} 1 {3} &> {4}'.format(
                masterdate.strftime('%Y%m%d')+'_crop',gc.rglks,gc.azlks,
                masterrslcdir,logfilename)
        rc = os.system(multicall)
        if rc != 0:
            print('Something went wrong multilooking the '\
                    'merged image. Log file {0}. Continuing with next '\
                    'acquisition date.'.format(logfilename), file=sys.stderr)
    #create file paths
        geocropdir = os.path.join(slaverslcdir,'geo')
        os.mkdir(geocropdir)
        demdir = os.path.join(procdir,'DEM')
        # rc = geocode_dem(masterrslcdir,geocropdir,demdir,
                # procdir,masterdate.strftime('%Y%m%d')+'_crop',gc.outres)
        
####
#####
        
############################################################ Coregester slave to master. See no missing bursts code
        logfile = os.path.join(procdir,'log','rdc_trans_'+
                               masterdate.strftime('%Y%m%d')+'_'+
                               slavedate.strftime('%Y%m%d')+'.log')
        print('Resampling cropped SLC using the existing LUT')
        if not SLC_interp_lt_S1_TOPS(slaveslctab,slavepar,croptab,
                                     cropfilename+'.par',lut,mastermlipar,
                                     slavemlipar,
                                     offfile,slaverslctab,slaverfilename,
                                     slaverfilename+'.par',logfile):
            print("\nError:", file=sys.stderr)
            print("Something went wrong resampling the slave SLC", file=sys.stderr)
            return 3

        logfilename =  os.path.join(procdir,'log',
                                            'multilookRSLC_{0}.log'.format(
                                                slavedate.strftime('%Y%m%d')))
        
    
############################################################ Pad cropped data
        azoffsetcall = ['master2mastercrop_offset.sh',
                        os.path.join(masterrslcdir,
                                     masterdate.strftime('%Y%m%d')),
                        'rslc',
                        os.path.join(masterrslcdir,
                                     masterdate.strftime('%Y%m%d'))+'_crop',
                        'rslc']
        azoffset = subp.check_output(azoffsetcall).strip()

        rslccall = ['rslc2rslc.sh',
                    os.path.join(masterrslcdir,
                                     masterdate.strftime('%Y%m%d')),
                    os.path.join(slaverslcdir,
                                     slavedate.strftime('%Y%m%d')),
                    azoffset]
        rc = subp.check_call(rslccall)
        if rc > 0:
            print("\nError:", file=sys.stderr)
            print("Something went wrong zero padding cropped RSLC.", file=sys.stderr)
            return 4

        logfilename =  os.path.join(procdir,'log',
                                            'multilookRSLC_{0}_crop.log'.format(
                                                slavedate.strftime('%Y%m%d')))

    #multilook cropped padded data
        multicall = 'multilookRSLC {0} {1} {2} 1 {3} &> {4}'.format(
                slavedate.strftime('%Y%m%d'),gc.rglks,gc.azlks,slaverslcdir,logfilename)
        rc = os.system(multicall)
        if rc != 0:
            print('Something went wrong multilooking the resampled slave image. Log file {0}. Continuing with next acquisition date.'.format(logfilename), file=sys.stderr)

############################################################ remove lock file
    os.remove(slaveLockFile)

    return 0

def regenerate_offfile(offfile, coreg_qual, mastermlipar):
    if os.path.exists(offfile):
        print('the off file already exists, stopping regeneration')
        return 1
    print('Regenerating off file from existing coreg_qual')
    rangepixels = grep1('range_samples:',mastermlipar.replace('slc.mli.par','slc.par')).split(':')[1].strip()
    width = grep1('range_samples:',mastermlipar).split(':')[1].strip()
    length = grep1('azimuth_lines:',mastermlipar).split(':')[1].strip()
    try:
        total_az_wrt_rdc = grep1('Total azimuth offset',coreg_qual).split(':')[1].strip().split()[0]
    except:
        print('Total azimuth offset not estimated in '+coreg_qual+'. Cancelling recoreg')
        return 1
    az_crosscorronly = grep1('intensity_matching:',coreg_qual).split(':')[1].strip().split()[0]
    azoff = float(total_az_wrt_rdc) - float(az_crosscorronly)
    
    with open(offfile, 'w') as f:
        f.write('Gamma Interferometric SAR Processor (ISP)\n')
        f.write('Interferogram and Image Offset Parameter File\n')
        f.write('initial_range_offset:                    0\n')
        f.write('initial_azimuth_offset:                  0\n')
        f.write('slc1_starting_range_pixel:               0\n')
        f.write('number_of_slc_range_pixels:          {0}\n'.format(rangepixels))
        f.write('offset_estimation_starting_range:        0\n')
        f.write('offset_estimation_ending_range:          0\n')
        f.write('offset_estimation_range_samples:        32\n')
        f.write('offset_estimation_range_spacing:         0\n')
        f.write('offset_estimation_starting_azimuth:      0\n')
        f.write('offset_estimation_ending_azimuth:        0\n')
        f.write('offset_estimation_azimuth_samples:      64\n')
        f.write('offset_estimation_azimuth_spacing:       0\n')
        f.write('offset_estimation_window_width:         64\n')
        f.write('offset_estimation_window_height:        64\n')
        f.write('offset_estimation_threshold:          0.10\n')
        f.write('range_offset_polynomial:         0.00000   0.0000e+00   0.0000e+00   0.0000e+00   0.0000e+00   0.0000e+00\n')
        f.write('azimuth_offset_polynomial:    {0} 0.0000e+00 0.0000e+00 0.0000e+00 0.0000e+00 0.0000e+00\n'.format(str(azoff)))
        f.write('slc1_starting_azimuth_line:               0\n')
        f.write('interferogram_azimuth_lines:           {0}\n'.format(length))
        f.write('interferogram_width:                   {0}\n'.format(width))
        f.write('first_nonzero_range_pixel:                0\n')
        f.write('number_of_nonzero_range_pixels:        {0}\n'.format(width))
        f.write('interferogram_range_looks:               20\n')
        f.write('interferogram_azimuth_looks:              4\n')
        f.write('interferogram_range_pixel_spacing:      46.591240   m\n')
        f.write('interferogram_azimuth_pixel_spacing:    55.839120   m\n')
        f.write('resampled_range_pixel_spacing:           0.000000   m\n')
        f.write('resampled_azimuth_pixel_spacing:         0.000000   m\n')
        f.write('resampled_starting_ground_range:         0.00000   m\n')
        f.write('resampled_pixels_per_line:               0\n')
        f.write('resampled_number_of_lines:               0\n\n')
    return 0


def recoreg_slave_old(slavedate,slcdir,rslcdir,masterdate,framename,procdir,lq):
    """ Recoregister and resample slave to master geometry. This assumes lookup
    tables have already been created.
    """
    print('\nRecoregistering slave {0}...'.format(slavedate.date()))
    
    #get/create slave slc directory paths
    slaveslcdir = os.path.join(slcdir,slavedate.strftime('%Y%m%d'))
    slaverslcdir = os.path.join(rslcdir,slavedate.strftime('%Y%m%d'))
    if not os.path.exists(slaverslcdir):
        return 6 # slave rslc dir doesn't exist!!

    #Create lock file for coregistration
    slaveLockFile = slaverslcdir+slavedate.strftime('/%Y%m%d.lock')
    open(slaveLockFile,'a').close() # OSX doesn't like os.mknod

    #master slc paths
    masterslcdir = os.path.join(slcdir,masterdate.strftime('%Y%m%d'))
    masterrslcdir = os.path.join(rslcdir,masterdate.strftime('%Y%m%d'))

    #get slave bursts
    slavebursts = lq.get_frame_bursts_on_date(framename,slavedate)
    #get master bursts
    masterbursts = lq.get_frame_bursts_on_date(framename,masterdate)
    #Get sorted list of swaths
    swathlist = [x[0].split('_')[1] for x in masterbursts]
    swathlist = set(swathlist)
    swathlist = list(swathlist)
    swathlist.sort()

    #Check for bursts missing in slaves but present in master
    missingbursts = []
    for mb in masterbursts:
        if not mb in slavebursts:
            missingbursts.append(mb)


############################################################ Create input tab files

    #master tab file
    masterslctab = os.path.join(procdir,'tab',masterdate.strftime('%Y%m%d')+'_tab')
    masterfilename = os.path.join(masterrslcdir,masterdate.strftime('%Y%m%d')+'.rslc')
    rc, msg = make_SLC_tab(masterslctab,masterfilename,swathlist)
    if rc > 0:
        print('Something went wrong creating the master tab file...')
        return 1
        
    #slave slc tab file
    slaveslctab = os.path.join(procdir,'tab',slavedate.strftime('%Y%m%d')+'_tab')
    slavefilename = os.path.join(slaveslcdir,slavedate.strftime('%Y%m%d')+'.slc')
    rc, msg = make_SLC_tab(slaveslctab,slavefilename,swathlist)
    if rc > 0:
        print('Something went wrong creating the slave tab file...')
        return 1

    #resampled slave slc tab file
    slaverslctab = os.path.join(procdir,'tab',slavedate.strftime('%Y%m%d')+'R_tab')
    slaverfilename = os.path.join(slaverslcdir,slavedate.strftime('%Y%m%d')+'.rslc')
    rc, msg = make_SLC_tab(slaverslctab,slaverfilename,swathlist)
    if rc > 0:
        print('Something went wrong creating the slave resampled tab file...')
        return 1

    pair = os.path.join(slaverslcdir,masterdate.strftime('%Y%m%d')+'_'+
                           slavedate.strftime('%Y%m%d'))
    #offset parameter file path
    offfile = pair+'.off'
    refine2file = offfile+'.refine2'
############################################################ Geomatric coregistration
    #if no missing bursts....
    #missingbursts = False
    if not missingbursts:
        print('debug place 3')
        print('All bursts available, no recropping of master necessary...')
    #master mli param file path
        mastermlipar = os.path.join(masterslcdir,
                                 masterdate.strftime('%Y%m%d')+'.slc.mli.par')
    #master param file path
        masterpar = os.path.join(masterslcdir,
                                 masterdate.strftime('%Y%m%d')+'.slc.par')
    #slave mli param file path
        slavemlipar = os.path.join(slaveslcdir,
                                slavedate.strftime('%Y%m%d')+'.slc.mli.par')
    #master param file path
        slavepar = os.path.join(slaveslcdir,
                                slavedate.strftime('%Y%m%d')+'.slc.par')
    #DEM height file
        demhgt = os.path.join(procdir,'geo',masterdate.strftime('%Y%m%d')+'.hgt')
    #lookup table path
        lut = os.path.join(slaverslcdir,
                           masterdate.strftime('%Y%m%d')+'_'+
                           slavedate.strftime('%Y%m%d')+'.slc.mli.lt')
        if not os.path.exists(lut):
            return 7 # Couldn't find previous look up table
        logfile = os.path.join(procdir,'log','rdc_trans_'+
                               masterdate.strftime('%Y%m%d')+'_'+
                               slavedate.strftime('%Y%m%d')+'.log')

        print('Resampling image (again)...')
        if not SLC_interp_lt_S1_TOPS(slaveslctab,slavepar,masterslctab,
                                     masterpar,lut,mastermlipar,slavemlipar,
                                     refine2file,slaverslctab,slaverfilename,
                                     slaverfilename+'.par',logfile):
            print("\nError:", file=sys.stderr)
            print("Something went wrong resampling the slave SLC", file=sys.stderr)
            return 3
        #Multilook (average) resampled slc
        logfilename =  os.path.join(procdir,'log',
                                            'multilookRSLC_{0}_crop.log'.format(
                                                slavedate.strftime('%Y%m%d')))
        multicall = 'multilookRSLC {0} {1} {2} 1 {3} &> {4}'.format(
                slavedate.strftime('%Y%m%d'),gc.rglks,gc.azlks,slaverslcdir,logfilename)
        rc = os.system(multicall)
        if rc != 0:
            print('Something went wrong multilooking the '\
                    'resampled slave image. Log file {0}. Continuing with '\
                    'next acquisition date.'.format(logfilename), file=sys.stderr)

    else:
############################################################ Crop master to fit with smaller slave
        print('Missing bursts at either or both ends of the scene, recropping '\
                'master, and auxiliary slave if necessary...')
    #loop through swaths
        for iw in ['IW1','IW2','IW3']:
        #get master and slave bursts
            iwthisburst_m = sorted([b[0] for b in masterbursts if iw in b[0]])
            iwthisburst_s = sorted([b[0] for b in slavebursts if iw in b[0]])
        #find indices which match slave bursts
            iwstart = iwthisburst_m.index(iwthisburst_s[0])+1
            iwstop = iwthisburst_m.index(iwthisburst_s[-1])+1

        #Create master slc tab file
            masterslctab = os.path.join(procdir,'tab',masterdate.strftime('%Y%m%d')
                                                                            +'_tab')
            masterslcfilename = os.path.join(masterrslcdir,masterdate.strftime('%Y%m%d')
                                                                            +'.rslc')
            rc, msg = make_SLC_tab(masterslctab,masterslcfilename,[iw])
            if rc > 0:
                print('Something went wrong creating the master tab file...')
                return 1
        #create a corresponding cropped master slc tab file
            croptab = os.path.join(procdir,'tab',masterdate.strftime('%Y%m%d')+
                                                                        'crop_tab')
            cropfilename = os.path.join(masterrslcdir,masterdate.strftime('%Y%m%d')
                                                                        +'_crop.rslc')
            rc, msg = make_SLC_tab(croptab,cropfilename,[iw])
            if rc > 0:
                print('Something went wrong creating the cropped master tab file...')
                return 1

            bursttab = os.path.join( procdir,'tab', 'master_crop_{0}_burst_tab'.format( iw ) )
            rc, msg = make_burst_tab( bursttab, iwstart, iwstop)
            if rc > 0:
                print('\nProblem creating burst '\
                                    'tab {0} for date {1}. Error '\
                                    'message:'\
                                    '\n {2}'\
                                    '\nContinuing with next '\
                                    'date.'.format( bursttab, 
                                            date, msg ), file=sys.stderr)
                os.remove(slaveLockFile)
                return 2
            logfilename = os.path.join(procdir,'log','SLC_copy_S1_TOPS_'+
                                       masterdate.strftime('%Y%m%d')+
                                       '_crop.log')
        #copy out bursts which are relevant to the cropped master
            if SLC_copy_S1_TOPS(masterslctab,croptab,bursttab,'',procdir,
                                                                    logfilename):
                print("\nError:", file=sys.stderr)
                print("Something went wrong recropping the master SLC", file=sys.stderr)
                os.remove(slaveLockFile)
                return 5

        logfilename = os.path.join(procdir,'log','SLC_mosaic_S1_TOPS_'+
                                   masterdate.strftime('%Y%m%d')+
                                   '_{0}_crop.log'.format(iw))

    #mosiac cropped file (does this mean if aux slave is present we only use that cropped file?)
        rc, msg = make_SLC_tab(croptab,cropfilename,swathlist)
        if not SLC_mosaic_S1_TOPS(croptab,cropfilename,gc.rglks,gc.azlks,logfilename):
            print('Something went wrong mosaicing subswaths '\
                    'together. Log file {0}. Continuing with next acquisition '\
                    'date.'.format(logfilename), file=sys.stderr)
            os.remove(slaveLockFile)
            return 3
        logfilename =  os.path.join(procdir,'log',
                            'multilookRSLC_{0}_crop.log'.format(
                                masterdate.strftime('%Y%m%d')))
    #multilook (average) cropped data
        multicall = 'multilookRSLC {0} {1} {2} 1 {3} &> {4}'.format(
                masterdate.strftime('%Y%m%d')+'_crop',gc.rglks,gc.azlks,
                masterrslcdir,logfilename)
        rc = os.system(multicall)
        if rc != 0:
            print('Something went wrong multilooking the '\
                    'merged image. Log file {0}. Continuing with next '\
                    'acquisition date.'.format(logfilename), file=sys.stderr)
    #create file paths
        geocropdir = os.path.join(slaverslcdir,'geo')
        os.mkdir(geocropdir)
        demdir = os.path.join(procdir,'DEM')
        # rc = geocode_dem(masterrslcdir,geocropdir,demdir,
                # procdir,masterdate.strftime('%Y%m%d')+'_crop',gc.outres)
        
####
#####
        
############################################################ Coregester slave to master. See no missing bursts code
        mastermlipar = os.path.join(masterrslcdir,
                                 masterdate.strftime('%Y%m%d')+'_crop.rslc.mli.par')
        masterpar = os.path.join(masterrslcdir,
                                 masterdate.strftime('%Y%m%d')+'_crop.slc.par')
        slavemlipar = os.path.join(slaveslcdir,
                                slavedate.strftime('%Y%m%d')+'.slc.mli.par')
        slavepar = os.path.join(slaveslcdir,
                                slavedate.strftime('%Y%m%d')+'.slc.par')
        demhgt = os.path.join(slaverslcdir,'geo',masterdate.strftime('%Y%m%d')
                                                                +'_crop.hgt')
        lut = os.path.join(slaverslcdir,
                           masterdate.strftime('%Y%m%d')+'_'+
                           slavedate.strftime('%Y%m%d')+'.slc.mli.lt')
        if not os.path.exists(lut):
            os.remove(slaveLockFile)
            return 7 # Couldn't find previous look up table
        logfile = os.path.join(procdir,'log','rdc_trans_'+
                               masterdate.strftime('%Y%m%d')+'_'+
                               slavedate.strftime('%Y%m%d')+'.log')
        print('Resampling image (again)...')
        if not SLC_interp_lt_S1_TOPS(slaveslctab,slavepar,croptab,
                                     cropfilename+'.par',lut,mastermlipar,
                                     slavemlipar,
                                     refine2file,slaverslctab,slaverfilename,
                                     slaverfilename+'.par',logfile):
            print("\nError:", file=sys.stderr)
            print("Something went wrong resampling the slave SLC", file=sys.stderr)
            os.remove(slaveLockFile)
            return 3

        logfilename =  os.path.join(procdir,'log',
                                            'multilookRSLC_{0}.log'.format(
                                                slavedate.strftime('%Y%m%d')))
        
    
############################################################ Pad cropped data
        azoffsetcall = ['master2mastercrop_offset.sh',
                        os.path.join(masterrslcdir,
                                     masterdate.strftime('%Y%m%d')),
                        'rslc',
                        os.path.join(masterrslcdir,
                                     masterdate.strftime('%Y%m%d'))+'_crop',
                        'rslc']
        azoffset = subp.check_output(azoffsetcall).strip()

        rslccall = ['rslc2rslc.sh',
                    os.path.join(masterrslcdir,
                                     masterdate.strftime('%Y%m%d')),
                    os.path.join(slaverslcdir,
                                     slavedate.strftime('%Y%m%d')),
                    azoffset]
        rc = subp.check_call(rslccall)
        if rc > 0:
            print("\nError:", file=sys.stderr)
            print("Something went wrong zero padding cropped RSLC.", file=sys.stderr)
            os.remove(slaveLockFile)
            return 4

        logfilename =  os.path.join(procdir,'log',
                                            'multilookRSLC_{0}_crop.log'.format(
                                                slavedate.strftime('%Y%m%d')))

    #multilook cropped padded data
        multicall = 'multilookRSLC {0} {1} {2} 1 {3} &> {4}'.format(
                slavedate.strftime('%Y%m%d'),gc.rglks,gc.azlks,slaverslcdir,logfilename)
        rc = os.system(multicall)
        if rc != 0:
            print('Something went wrong multilooking the resampled slave image. Log file {0}. Continuing with next acquisition date.'.format(logfilename), file=sys.stderr)

############################################################ remove lock file
    os.remove(slaveLockFile)

    return 0




################################################################################
# Co-register slave function for stripmaps
# (by Daniel Juncu, Uni of Leeds)
# currently co-registration follows ISP_users_guide A.3, p.47 ff
# alternative approach could be implemented from GEO_users_guide Section I, p.81 ff
################################################################################
def coreg_slave_common_sm(procdir,masterdate,masterrslcdir,slavedate,slaveslcdir,slaverslcdir,slaverfilename,qualityfile,crop = False):
    if not crop:
        croptext=''
        geodir = os.path.join(procdir,'geo')
    else:
        croptext='_crop_'+slavedate.strftime('%Y%m%d')
        geodir = os.path.join(slaverslcdir,'geo')
    mastermli = os.path.join(procdir,'SLC',masterdate.strftime('%Y%m%d'),
                             masterdate.strftime('%Y%m%d')+croptext+'.slc.mli')
    #master mli param file path
    masterslctab = os.path.join(procdir,'tab',masterdate.strftime('%Y%m%d')+croptext+'_tab')
    mastermlipar = os.path.join(procdir,'SLC',masterdate.strftime('%Y%m%d'),
                             masterdate.strftime('%Y%m%d')+croptext+'.slc.mli.par')
    [mliwidth,mlilength]=get_mli_size(mastermlipar)
    #master param file path
    masterpar = os.path.join(masterrslcdir,
                             masterdate.strftime('%Y%m%d')+croptext+'.rslc.par')
    masterslcfilename = os.path.join(masterrslcdir,masterdate.strftime('%Y%m%d')+croptext+'.rslc')
    slaveslcfilename    = os.path.join(slaveslcdir,slavedate.strftime('%Y%m%d')+croptext+'.slc')
    #slave mli param file path
    slaveRmlipar = os.path.join(slaverslcdir,
                            slavedate.strftime('%Y%m%d')+'.rslc.mli.par')
    slaveRmlifile = os.path.join(slaverslcdir,
                            slavedate.strftime('%Y%m%d')+'.rslc.mli')
    #master param file path
    slavepar = os.path.join(slaveslcdir,
                            slavedate.strftime('%Y%m%d')+'.slc.par')
    slaveRpar = os.path.join(slaverslcdir,
                            slavedate.strftime('%Y%m%d')+'.rslc.par')
    slaverslctab = os.path.join(procdir,'tab',slavedate.strftime('%Y%m%d')+'R_tab')
    #DEM height file
    demhgt = os.path.join(geodir,masterdate.strftime('%Y%m%d')+croptext+'.hgt')
    #lookup table path
    lut = os.path.join(slaverslcdir,
                       masterdate.strftime('%Y%m%d')+'_'+
                       slavedate.strftime('%Y%m%d')+'.slc.mli.lt')
    pair = os.path.join(slaverslcdir,masterdate.strftime('%Y%m%d')+'_'+
                           slavedate.strftime('%Y%m%d'))
    #offset parameter file path
    offfile = pair+'.off'
    dofffile = pair+'.doff'
    #setting logfile and going on
    logfile_offset = os.path.join(procdir,'log','create_offset_'+
                           masterdate.strftime('%Y%m%d')+'_'+
                           slavedate.strftime('%Y%m%d')+'.log')
    if not create_offset(masterpar,slavepar,offfile,str(1),
                         str(1),logfile_offset):
        print("\nError:", file=sys.stderr)
        print("Something went wrong creating the offset file.", file=sys.stderr)
        return 4
    with open(qualityfile, "a") as myfile:
            myfile.write("Improvement of coregistration offset using intensity cross-correlation:\n")
    print("Offset estimation using intensity cross-correlation")
    shutil.copyfile(offfile,offfile+'.start')
    print("Offset estimation using using orbit information")
    # estimation of offset using orbit information
    if not init_offset_orbit(masterpar,slavepar,offfile,logfile_offset):
        print("\nError:", file=sys.stderr)
        print("Something went wrong estimating offset using orbit information (1).", file=sys.stderr)
        return 4
    print("Improve offset estimate")
    # improve estimate (multi-looked)
    #previously: gc.rglks, gc. azlks
    ##if not init_offset(masterslcfilename,slaveslcfilename,masterpar,slavepar,offfile, 4, 4, logfile_offset, 512, 512, 0.15, False):
    ##    print("\nError:", file=sys.stderr)
    ##    print("Something went wrong estimating offset using orbit information (2).", file=sys.stderr)
    ##    return 4
    ### improve estimate (single-look)
    ##if not init_offset(masterslcfilename,slaveslcfilename,masterpar,slavepar,offfile, 1, 1, logfile_offset, 128, 128, 0.15, False):
    ##    print("\nError:", file=sys.stderr)
    ##    print("Something went wrong estimating offset using orbit information (3).", file=sys.stderr)
    ##    return 4
    logfile = os.path.join(procdir,'log','offset_pwr_'+
                       masterdate.strftime('%Y%m%d')+'_'+
                       slavedate.strftime('%Y%m%d')+'.log')
    if not offset_pwr(masterslcfilename,slaveslcfilename,masterpar,slavepar,offfile,pair,logfile):
        print("\nError:", file=sys.stderr)
        print("Something went wrong with offset estimation.", file=sys.stderr)
        return 4
    #then we do offset fitting (gamma recommends offset_fit and not offset_fitm)
    if not offset_fit(pair,offfile,pair+'.off.out',acqMode='sm'):
        print("\nError:", file=sys.stderr)
        print("Something went wrong with offset fitting.", file=sys.stderr)
        return 4
    fittmp = grep1('final model fit std. dev.',pair+'.off.out')
    range_stdev = float(fittmp.split(':')[1].split()[0])
    azimuth_stdev = float(fittmp.split(':')[2])
    daztmp = grep1('azimuth_offset_polynomial:',offfile)
    drtmp = grep1('range_offset_polynomial:',offfile)
    daz = float(daztmp.split()[1])
    daz_mli = daz/gc.azlks
    dr = float(drtmp.split()[1])
    dr_mli = dr/gc.rglks
    if os.path.exists(pair+'.diff_par'): os.remove(pair+'.diff_par')
    logfile = os.path.join(procdir,'log',
                   'resampling_'+
                   masterdate.strftime('%Y%m%d')+'_'+
                   slavedate.strftime('%Y%m%d')+'.log')
    # resample 
    print("Resampling...")
    if not SLC_interp(slaveslcfilename, masterpar, slavepar, offfile, slaverfilename, slaveRpar,logfile):
        print("\nError:", file=sys.stderr)
        print("Something went wrong resampling the slave SLC", file=sys.stderr)
        return 3
    # multi look RSLC
    logfile = os.path.join(procdir,'log',
                   'multilooking'+
                   masterdate.strftime('%Y%m%d')+'_'+
                   slavedate.strftime('%Y%m%d')+'.log')
    print("multi look RSLC...")    
    if not multi_look(slaverfilename, slaveRpar, slaveRmlifile, slaveRmlipar, gc.rglks, gc.azlks, logfile):
        print("\nError:", file=sys.stderr)
        print("Something went wrong multilooking", file=sys.stderr)
        return 3
    logfile = os.path.join(procdir,'log','rdc_trans_'+
                           masterdate.strftime('%Y%m%d')+'_'+
                           slavedate.strftime('%Y%m%d')+'.log')
    # create lookup table from master and slave mli par and master height file
    print('Derive lookup table for MLI coregistration...')
    if not rdc_trans(mastermlipar,demhgt,slaveRmlipar,lut,logfile):
        print("\nError:", file=sys.stderr)
        print("Something went wrong during lookup table generation.", file=sys.stderr)
        return 2    
    return 0
#########################################################################

#########################################################################
# coreg slave to master for stripmaps
####
def coreg_slave_sm(slavedate,slcdir,rslcdir,masterdate,framename,procdir, lq, job_id):
    """ Coregister and resample slave to master geometry
    """ 
    print('\nCoregistering slave {0}...'.format(slavedate.date()))
    #get/create slave slc directory paths
    slaveslcdir = os.path.join(slcdir,slavedate.strftime('%Y%m%d'))
    slaverslcdir = os.path.join(rslcdir,slavedate.strftime('%Y%m%d'))
    if os.path.exists(slaverslcdir):
        shutil.rmtree(slaverslcdir)
    os.mkdir(slaverslcdir)
    #Create lock file for coregistration
    slaveLockFile = slaverslcdir+slavedate.strftime('/%Y%m%d.lock')
    open(slaveLockFile,'a').close() # OSX doesn't like os.mknod
    #master slc paths
    masterslcdir = os.path.join(slcdir,masterdate.strftime('%Y%m%d'))
    masterrslcdir = os.path.join(rslcdir,masterdate.strftime('%Y%m%d'))
    #calc diff date
    masterdiff = abs(slavedate.date()-masterdate)
    print("master difference is {0}".format(masterdiff))
    #create slave diff list
    procslavelist = [s for s in os.listdir(rslcdir) if 
            s != masterdate.strftime( '%Y%m%d') and s != slavedate.strftime('%Y%m%d')]
    print("Found slaves ... {0}".format(procslavelist))
    procslavediff = [dt.datetime.strptime(s,'%Y%m%d')-slavedate for s in
                                                                    procslavelist]
############################################################ Create input tab files
    masterfilename = os.path.join(masterrslcdir,masterdate.strftime('%Y%m%d')+'.rslc')
    slavefilename = os.path.join(slaveslcdir,slavedate.strftime('%Y%m%d')+'.slc')
    slaverfilename = os.path.join(slaverslcdir,slavedate.strftime('%Y%m%d')+'.rslc')
############################################################ Geometric coregistration
    qualityfile=os.path.join(procdir,'log','coreg_quality_'+
                              masterdate.strftime('%Y%m%d')+'_'+
                              slavedate.strftime('%Y%m%d')+'.log')
    rc = coreg_slave_common_sm(procdir,masterdate,masterrslcdir,slavedate,slaveslcdir,slaverslcdir,slaverfilename,qualityfile)
    if rc != 0:
        print("\nError:", file=sys.stderr)
        print("Something went wrong during coregistration", file=sys.stderr)
        return rc
############################################################ remove lock file
    os.remove(slaveLockFile)
############################################################ Update job database (if not in batch mode)
    # Populate the RSLC table with the location of the new slave RSLC
    print("RSLC succesfully made, now to write to RSLC table...")
    print("job_id: %d" % job_id)
    if job_id != -1:
        print('Attempting to input new rslc SLAVE product')
        lq.set_new_rslc_product(job_id,slavedate.strftime('%Y%m%d'),
                                masterrslcdir+'/'+masterdate.strftime('%Y%m%d')+'.rslc',
                                slavedate.strftime('%Y%m%d') + '.rslc', slaverslcdir)
    return 0
