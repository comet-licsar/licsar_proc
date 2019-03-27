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

from glob import glob

#Import Gamma functions
from gamma_functions import *
#And import the global configuration
import global_config as gc

################################################################################
#Link master rslc function
################################################################################
def link_master_rslc(masterslcdir, rslcdir, masterdate, lq, job_id):
    """ Creates RSLC dir and links master SLC images
    """

############################################################Create directories
    if not os.path.exists(rslcdir):
        os.mkdir(rslcdir)
    
    masterrslcdir = os.path.join(rslcdir,masterdate.strftime('%Y%m%d'))
    if os.path.exists(masterrslcdir):
        shutil.rmtree(masterrslcdir)
    os.mkdir(masterrslcdir)
    
############################################################ Iterate through 
                                                        #and syslink master slcs
    filelist = [x for x in os.listdir(masterslcdir) if '{0}.'.format(
                                    masterdate.strftime('%Y%m%d')) in x]
    filelistnew = [x.replace('.slc','.rslc') for x in filelist]
    try:
        for f1,f2 in zip(filelist,filelistnew):
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
def geocode_dem(masterslcdir,geodir,demdir,procdir,masterdate,outres):
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
    mli2 = mlipar[:-4] #The original master mli
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
def coreg_slave(slavedate,slcdir,rslcdir,masterdate,framename,procdir, lq, job_id):
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
    # print "proc slave diffs {0}".format(procslavediff)
    #get slave bursts
    slavebursts = lq.get_frame_bursts_on_date(framename,slavedate)
    #get master bursts
    masterbursts = lq.get_frame_bursts_on_date(framename,masterdate)

############################################################ Find nearest (date) coregestered slave (already processed slave) to use as auxiliary
    cond = True
    while cond: #sweep through processed slaves till an appropriate aux is found, or till no processed slave available
        if procslavediff: # other slaves have been processed
            if np.min(np.abs(procslavediff)) < masterdiff and np.min(
                    np.abs(procslavediff)) > dt.timedelta(0): 
                # slaves closer than master
                slave3ix = np.argsort(np.abs(procslavediff))[0]
                slave3date = dt.datetime.strptime(procslavelist[slave3ix],'%Y%m%d')
                print("found potential aux slave date {0}".format(slave3date))
                slave3rslcdir = os.path.join(rslcdir,slave3date.strftime('%Y%m%d'))
                print(slavebursts)
                print(masterbursts)
                slave3bursts = masterbursts
                #slave3bursts = lq.get_frame_bursts_on_date(framename,slave3date)
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
                    procslavelist.pop(slave3ix)
                    procslavediff.pop(slave3ix)
                    print("not valid as aux slave")
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

############################################################ Geomatric coregistration
    #if no missing bursts....
    if not missingbursts:
        print('All bursts available, no recropping of master necessary...')
        mastermli = os.path.join(masterslcdir,
                                 masterdate.strftime('%Y%m%d')+'.slc.mli')
        #master mli param file path
        mastermlipar = os.path.join(masterslcdir,
                                 masterdate.strftime('%Y%m%d')+'.slc.mli.par')
        [mliwidth,mlilength]=get_mli_size(mastermlipar)
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
        logfile = os.path.join(procdir,'log','rdc_trans_'+
                               masterdate.strftime('%Y%m%d')+'_'+
                               slavedate.strftime('%Y%m%d')+'.log')

        #create lookup table from master and slave mli par and master height file
        print('Performing the geometric coregistration...')
        #print(mastermlipar,demhgt,slavemlipar,lut,logfile)
        #a=rdc_trans(mastermlipar,demhgt,slavemlipar,lut,logfile)
        if not rdc_trans(mastermlipar,demhgt,slavemlipar,lut,logfile):
            print("\nError:", file=sys.stderr)
            print("Something went wrong during geometric coregistration.", file=sys.stderr)
            return 2
        #logfile = os.path.join(procdir,'log','SLC_interp_lt_S1_TOPS_'+
        #                       masterdate.strftime('%Y%m%d')+'_'+
        #                       slavedate.strftime('%Y%m%d')+'.log')
        #interpolate slave image to master using lookup table

        ###########################################################################################################
        #ML: to reflect gamma 'recipe' from 20181130 (and previous but.. still not reflected version)
        # we will iterate already the mli-based improvements to reach 0.01 px
        pair = os.path.join(slaverslcdir,masterdate.strftime('%Y%m%d')+'_'+
                               slavedate.strftime('%Y%m%d'))
        #offset parameter file path
        offfile = pair+'.off'
        dofffile = pair+'.doff'
        #correction of offset estimation samples (max 64x64)
        [width,length]=get_mli_size(masterpar)
        rs=round(int(width)/64)
        az=round(int(length)/64)
        rstep = 64 if 64 > rs else rs
        azstep= 32 if 32 > az else az
        #setting logfile and going on
        logfile_offset = os.path.join(procdir,'log','create_offset_'+
                               masterdate.strftime('%Y%m%d')+'_'+
                               slavedate.strftime('%Y%m%d')+'.log')
        if not create_offset(masterpar,slavepar,offfile,str(gc.rglks),
                             str(gc.azlks),logfile_offset):
            print("\nError:", file=sys.stderr)
            print("Something went wrong creating the offset file.", file=sys.stderr)
            return 4
        qualityfile=os.path.join(procdir,'log','coreg_quality_'+
                              masterdate.strftime('%Y%m%d')+'_'+
                              slavedate.strftime('%Y%m%d')+'.log')
        with open(qualityfile, "a") as myfile:
                myfile.write("Iterative improvement of refinement offset using matching:")
        daz10000=10000
        it=0
        itmax=5
        while (daz10000 > 100 or daz10000 < -100) and (it < itmax):
            it += 1
            print("Offset refinement - iteration "+str(it))
            shutil.copyfile(offfile,offfile+'.start')
            print('Resampling image...')
            if not SLC_interp_lt_S1_TOPS(slaveslctab,slavepar,masterslctab,
                                         masterpar,lut,mastermlipar,slavemlipar,
                                         offfile+'.start',slaverslctab,slaverfilename,
                                         slaverfilename+'.par',logfile):
                print("\nError:", file=sys.stderr)
                print("Something went wrong resampling the slave SLC", file=sys.stderr)
                return 3
            #we are using some temporary dofffile here, not sure why but was in gamma..
            if os.path.exists(dofffile): os.remove(dofffile)
            if not create_offset(masterpar,slaverfilename+'.par',dofffile,str(gc.rglks),
                             str(gc.azlks),logfile_offset):
                print("\nError:", file=sys.stderr)
                print("Something went wrong creating the offset file.", file=sys.stderr)
                return 4
            #we do offset tracking between SLC images using intensity cross-correlation
            logfile = os.path.join(procdir,'log','offset_pwr_tracking_'+
                               masterdate.strftime('%Y%m%d')+'_'+
                               slavedate.strftime('%Y%m%d')+'.log')
            if not offset_pwr_tracking(masterfilename,slaverfilename,masterpar,slaverfilename+'.par',
                             dofffile,pair,rstep,azstep,logfile):
                print("\nError:", file=sys.stderr)
                print("Something went wrong with offset tracking.", file=sys.stderr)
                return 4
            #then we do offset fitting (gamma recommends offset_fit and not offset_fitm)
            if not offset_fit(pair,dofffile,pair+'.off.out.'+str(it)):
                print("\nError:", file=sys.stderr)
                print("Something went wrong with offset fitting.", file=sys.stderr)
                return 4
            fittmp = grep1('final model fit std. dev.',pair+'.off.out.'+str(it))
            range_stdev = float(fittmp.split(':')[1].split()[0])
            azimuth_stdev = float(fittmp.split(':')[2])
            daztmp = grep1('azimuth_offset_polynomial:',dofffile)
            drtmp = grep1('range_offset_polynomial:',dofffile)
            daz10000 = int(float(daztmp.split()[1])*10000)
            daz = float(daztmp.split()[1])
            daz_mli = daz/gc.azlks
            dr = float(drtmp.split()[1])
            dr_mli = dr/gc.rglks
            with open(pair+'.refinement.iteration.'+str(it), "w") as myfile:
                myfile.write("dr_mli: "+str(dr_mli)+"    daz_mli: "+str(daz_mli))
            if os.path.exists(pair+'.diff_par'): os.remove(pair+'.diff_par')
            logfile = os.path.join(procdir,'log',
                           'create_diff_par_'+
                           masterdate.strftime('%Y%m%d')+'_'+
                           slavedate.strftime('%Y%m%d')+'.log')
            if not create_diff_par(mastermlipar,mastermlipar,pair+'.diff_par','1','0',logfile):
                print("\nError:", file=sys.stderr)
                print("Something went wrong during the DIFF parameter file creation", file=sys.stderr)
                return 3
            shutil.copyfile(pair+'.diff_par',pair+'.diff_par.'+str(it))
            if not set_value(pair+'.diff_par',pair+'.diff_par',"range_offset_polynomial",
                           str(dr_mli)+"   0.0000e+00   0.0000e+00   0.0000e+00   0.0000e+00   0.0000e+00"):
                print("\nError:", file=sys.stderr)
                print("Something went wrong during the range value setting", file=sys.stderr)
                return 3
            if not set_value(pair+'.diff_par',pair+'.diff_par',"azimuth_offset_polynomial",
                           str(daz_mli)+"   0.0000e+00   0.0000e+00   0.0000e+00   0.0000e+00   0.0000e+00"):
                print("\nError:", file=sys.stderr)
                print("Something went wrong during the azimuth value setting", file=sys.stderr)
                return 3
            shutil.move(lut,lut+'.tmp.'+str(it))
            logfile = os.path.join(procdir,'log',
                           'gc_map_fine_{0}.log'.format(masterdate.strftime('%Y%m%d')))
            if not gc_map_fine(lut+'.tmp.'+str(it),str(mliwidth),pair+'.diff_par',lut,'1',logfile):
                # Error in lookup table refining
                print("\nError:", file=sys.stderr)
                print("Something went wrong refining the lookup table.", file=sys.stderr)
                return 4
            with open(qualityfile, "a") as myfile:
                myfile.write("matching_iteration_"+str(it)+": "+str(daz)+" "+str(dr)+" "+str(daz_mli)+" "+str(dr_mli)+" (daz dr  daz_mli dr_mli)")
                myfile.write("matching_iteration_"+str(it)+": "+str(azimuth_stdev)+" "+str(range_stdev) +"  (azi/rg std dev) ")
        #
        # Iterative improvement of azimuth refinement using spectral diversity method
        #
        #
        ####### These lines may be used just once? Only for master, before all the coregs start
        mazpoly=os.path.join(slaverslcdir,masterdate.strftime('%Y%m%d')+'.az_ovr.poly')
        if os.path.exists(mazpoly): os.remove(mazpoly)
        logfile = os.path.join(procdir,'log',
                       'S1_poly_overlap_{0}.log'.format(masterdate.strftime('%Y%m%d')
                       +'_'+slavedate.strftime('%Y%m%d')))
        if not S1_poly_overlap(masterslctab,str(gc.rglks),str(gc.azlks),mazpoly,'1',logfile):
            # Error in lookup table refining
            print("\nError:", file=sys.stderr)
            print("Something went wrong using S1_poly_overlap", file=sys.stderr)
            return 4
        #
        #[mastermliwidth,mastermlilength]=get_mli_size(mastermlipar)
        mazovr=os.path.join(slaverslcdir,masterdate.strftime('%Y%m%d')+'.az_ovr')
        logfile = os.path.join(procdir,'log',
                       'polymath_{0}.log'.format(masterdate.strftime('%Y%m%d')
                       +'_'+slavedate.strftime('%Y%m%d')))
        if not poly_math(mastermli,mazovr,str(mliwidth),mazpoly,logfile):
            # Error in lookup table refining
            print("\nError:", file=sys.stderr)
            print("Something went wrong using S1_poly_overlap", file=sys.stderr)
            return 4
        # maybe should be here?
        logfile = os.path.join(procdir,'log',
                       'raspwr_{0}.log'.format(masterdate.strftime('%Y%m%d')
                       +'_'+slavedate.strftime('%Y%m%d')))
        if not raspwr(mazovr,str(mliwidth),'1','1',mazovr+'.ras',logfile):
            print("\nError:", file=sys.stderr)
            print("Something went wrong using raspwr", file=sys.stderr)
            return 4
        lutazovr=os.path.join(slaverslcdir,masterdate.strftime('%Y%m%d')+'.lt.az_ovr')
        logfile = os.path.join(procdir,'log',
                       'mask_class_{0}.log'.format(masterdate.strftime('%Y%m%d')
                       +'_'+slavedate.strftime('%Y%m%d')))
        if not mask_class(mazovr+'.ras',lut,lutazovr,logfile):
            print("\nError:", file=sys.stderr)
            print("Something went wrong using mask_class", file=sys.stderr)
            return 4
        ####### end of what should be done only for master
        #
        #
        #
        daz10000=10000
        it=0
        itmax=5
        with open(qualityfile, "a") as myfile:
                myfile.write("Iterative improvement of refinement offset azimuth overlap regions:")
        # iterate while azimuth correction >= 0.0005 SLC pixel
        while (daz10000 > 5 or daz10000 < -5) and (it < itmax):
            it += 1
            print('offset refinement using spectral diversity in azimuth overlap region iteration '+str(it))
            shutil.copyfile(offfile,offfile+'.start')
            logfile_interp_lt2 = os.path.join(procdir,'log',
                       'SLC_interp_lt_ScanSAR.{0}.2.out'.format(slavedate.strftime('%Y%m%d')))
            if not SLC_interp_lt_S1_TOPS(slaveslctab,slavepar,masterslctab,
                                         masterpar,lutazovr,mastermlipar,slavemlipar,
                                         offfile+'.start',slaverslctab,slaverfilename,
                                         slaverfilename+'.par',logfile_interp_lt2):
                print("\nError:", file=sys.stderr)
                print("Something went wrong resampling the slave SLC", file=sys.stderr)
                return 3
            offazovrout=offfile+'.az_ovr.'+str(it)+'.out'
            print('Getting offsets using spectral diversity...')
            logfile = os.path.join(procdir,'log','S1_coreg_overlap_'+
                                       masterdate.strftime('%Y%m%d')+'_'+
                                       slavedate.strftime('%Y%m%d')+'.'+str(it)+'.log')
            #note that the slave3tab is already set - and is empty string if should not be used
            if not S1_coreg_overlap(masterslctab,slaverslctab,masterdate.strftime('%Y%m%d')+'_'+slavedate.strftime('%Y%m%d'),
                                                offfile+'.start',offfile,slave3tab,offazovrout):
                print("\nError:", file=sys.stderr)
                print("Something went wrong during the first offset refinement using spectral diversity.", file=sys.stderr)
                return 4
            with open(qualityfile, "a") as myfile:
                myfile.write("az_ovr_iteration_"+str(it)+": "+str(daz)+" (daz in SLC pixel)")
            daztmp = grep1('azimuth_pixel_offset',offazovrout)
            daz = float(daztmp.split()[1])      
            daz10000 = int(daz*10000)
            shutil.copyfile(offfile,offfile+'.az_ovr.'+str(it))
        #
        #resample full data set
        print('Resampling the full data set')
        logfile_interp_lt3 = os.path.join(procdir,'log',
                       'SLC_interp_lt_ScanSAR.{0}.3.out'.format(slavedate.strftime('%Y%m%d')))
        if not SLC_interp_lt_S1_TOPS(slaveslctab,slavepar,masterslctab,
                                     masterpar,lut,mastermlipar,slavemlipar,
                                     offfile,slaverslctab,slaverfilename,
                                     slaverfilename+'.par',logfile_interp_lt3):
            print("\nError:", file=sys.stderr)
            print("Something went wrong resampling the slave SLC", file=sys.stderr)
            return 3
#########################################################################
        #Do the RSLC mosaic
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

    else:
############################################################ Crop master to fit with smaller slave
        # this solution was working also with the 20181130 gamma codes. Keeping as it is..
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
                rc, msg = make_SLC_tab(slave3slctab,slave3slcfilename,[iw])
                if rc > 0:
                    print('Something went wrong creating the auxiliary slave tab file...')
                    return 1
            #create a corresponding crop tab file for aux slave
                slave3_croptab = os.path.join(procdir,'tab',
                                slave3date.strftime('%Y%m%d')+'crop_tab')
                slave3_cropfilename = os.path.join(slave3rslcdir,
                                slave3date.strftime('%Y%m%d')+'_crop.rslc')
                rc, msg = make_SLC_tab(slave3_croptab,slave3_cropfilename,[iw])
                if rc > 0:
                    print('Something went wrong creating the cropped '\
                                        'auxiliary slave tab file...')
                    return 1
                logfilename = os.path.join(procdir,'log','SLC_copy_S1_TOPS_'+
                                           slave3date.strftime('%Y%m%d')+
                                           '_crop.log')
                bursttab = os.path.join( procdir,'tab', 'master_aux_crop_{0}_burst_tab'.format( iw ) )
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
                if SLC_copy_S1_TOPS(masterslctab,slave3_croptab,
                                   bursttab,'',procdir,logfilename):
                    print("\nError:", file=sys.stderr)
                    print("Something went wrong recropping the master SLC", file=sys.stderr)
                    return 5
                    
        #mosaic cropped files
        logfilename = os.path.join(procdir,'log','SLC_mosaic_S1_TOPS_'+
                                   masterdate.strftime('%Y%m%d')+
                                   '_{0}_crop.log'.format(iw))
        rc, msg = make_SLC_tab(croptab,cropfilename,swathlist)
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
                            'multilookRSLC_{0}_crop.log'.format(
                                masterdate.strftime('%Y%m%d')))

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
        rc = geocode_dem(masterrslcdir,geocropdir,demdir,
                procdir,masterdate.strftime('%Y%m%d')+'_crop',gc.outres)
        
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
        logfile = os.path.join(procdir,'log','rdc_trans_'+
                               masterdate.strftime('%Y%m%d')+'_'+
                               slavedate.strftime('%Y%m%d')+'.log')
        print('\nPerforming the geometric coregistration...')
        if not rdc_trans(mastermlipar,demhgt,slavemlipar,lut,logfile):
            print("\nError:", file=sys.stderr)
            print("Something went wrong during geometric coregistration.", file=sys.stderr)
            return 2
        logfile = os.path.join(procdir,'log','SLC_interp_lt_S1_TOPS_'+
                               masterdate.strftime('%Y%m%d')+'_'+
                               slavedate.strftime('%Y%m%d')+'.log')
        print('Resampling image...')
        if not SLC_interp_lt_S1_TOPS(slaveslctab,slavepar,croptab,
                                     cropfilename+'.par',lut,mastermlipar,
                                     slavemlipar,'-',slaverslctab,
                                     slaverfilename,slaverfilename+'.par',
                                     logfile):
            print("\nError:", file=sys.stderr)
            print("Something went wrong resampling the slave SLC", file=sys.stderr)
            return 3
        pair = os.path.join(slaverslcdir,masterdate.strftime('%Y%m%d')+'_'+
                               slavedate.strftime('%Y%m%d'))
        offfile = pair+'.off'
        logfile = os.path.join(procdir,'log','create_offset_'+
                               masterdate.strftime('%Y%m%d')+'_'+
                               slavedate.strftime('%Y%m%d')+'.log')
        if not create_offset(cropfilename+'.par',slaverfilename+'.par',
                             offfile,str(gc.rglks),
                             str(gc.azlks),logfile):
            print("\nError:", file=sys.stderr)
            print("Something went wrong creating the offset file.", file=sys.stderr)
            return 4
        refine1file = offfile+'.refine1'
        logfile = os.path.join(procdir,'log','S1_coreg_overlap_'+
                               masterdate.strftime('%Y%m%d')+'_'+
                               slavedate.strftime('%Y%m%d')+'-1.log')
        print('Getting offsets using spectral diversity...')
#        if not S1_coreg_overlap(slave3tab,slaverslctab,pair,
#        print 'S1_coreg_overlap(' + croptab + ',' + slaverslctab + ',' + pair + ',' + offfile + ',' + refine1file + ',' + slave3tab + ',' + logfile + ')'
    #This is to correct possibility of having cropped master for coregistration and no aux slave
    #    slave3tab_crop = ''
    #if slave3tab:
    #        slave3tab_crop = croptab
        if not S1_coreg_overlap(croptab,slaverslctab,pair,
                        offfile,refine1file,slave3_croptab,logfile):
#                        offfile,refine1file,slave3tab,logfile):
            print("\nError:", file=sys.stderr)
            print("Something went wrong during the first "\
                    "offset refinement using spectral diversity.", file=sys.stderr)
            return 4      

        logfile = os.path.join(procdir,'log','SLC_interp_lt_S1_TOPS_'+
                               masterdate.strftime('%Y%m%d')+'_'+
                               slavedate.strftime('%Y%m%d')+'_refine1.log')
        print('Resampling image...')
        if not SLC_interp_lt_S1_TOPS(slaveslctab,slavepar,croptab,
                                     cropfilename+'.par',lut,mastermlipar,
                                     slavemlipar,
                                     refine1file,slaverslctab,slaverfilename,
                                     slaverfilename+'.par',logfile):
            print("\nError:", file=sys.stderr)
            print("Something went wrong resampling the slave SLC", file=sys.stderr)
            return 3

        refine2file = offfile+'.refine2'
        logfile = os.path.join(procdir,'log','S1_coreg_overlap_'+
                               masterdate.strftime('%Y%m%d')+'_'+
                               slavedate.strftime('%Y%m%d')+'-2.log')
        print('Getting offsets using spectral diversity...')
#        if not S1_coreg_overlap(slave3tab,slaverslctab,pair,refine1file,
        if not S1_coreg_overlap(croptab,slaverslctab,pair,refine1file,
                refine2file,slave3_croptab,logfile):
#                refine2file,slave3tab,logfile):
            print("\nError:", file=sys.stderr)
            print("Something went wrong during the first "\
                    "offset refinement using spectral diversity.", file=sys.stderr)
            return 4      

        logfile = os.path.join(procdir,'log','SLC_interp_lt_S1_TOPS_'+
                               masterdate.strftime('%Y%m%d')+'_'+
                               slavedate.strftime('%Y%m%d')+'_refine2.log')
        print('Resampling image (again)...')
        if not SLC_interp_lt_S1_TOPS(slaveslctab,slavepar,croptab,
                                     cropfilename+'.par',lut,mastermlipar,
                                     slavemlipar,
                                     refine2file,slaverslctab,slaverfilename,
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
    if not missingbursts:
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

    #mosiac cropped file (does this mean if aux slave is present we only use that cropped file?)
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
