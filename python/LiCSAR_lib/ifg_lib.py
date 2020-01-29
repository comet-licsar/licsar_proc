################################################################################
#Imports
################################################################################
import os
import shutil
import sys
import datetime as dt
import glob

from LiCSAR_misc import is_non_zero_file, grep1
from LiCSAR_lib.coreg_lib import get_mli_size
from gamma_functions import *
import global_config as gc

################################################################################
#Make interferograms functions
################################################################################
def make_interferogram(origmasterdate,masterdate,slavedate,procdir, lq, job_id, rglks = gc.rglks, azlks = gc.azlks):
    """
    Makes a singular interferogram between the SLC on given master and slave dates.
    """
############################################################ Get orig. master slc
    mastermli = os.path.join(procdir,'SLC',
                             origmasterdate.strftime('%Y%m%d'),
                             origmasterdate.strftime('%Y%m%d')+
                             '.slc.mli.par')
    mastertab = os.path.join(procdir,'tab',
                             origmasterdate.strftime('%Y%m%d')
                             +'_tab')
    #read width and length
    [width, length] = get_mli_size(mastermli)

############################################################ Create interferogram name
    pair = '{0}_{1}'.format(masterdate.strftime('%Y%m%d'),
                            slavedate.strftime('%Y%m%d'))
    print('Computing interferogram {0}'.format(pair))

############################################################ Create a log with basic ifg info
    qualityfile=os.path.join(procdir,'log','ifg_quality_'+
                              masterdate.strftime('%Y%m%d')+'_'+
                              slavedate.strftime('%Y%m%d')+'.log')
    with open(qualityfile, "a") as myfile:
            myfile.write("Basic interferogram information: \n")
    
############################################################ Create directory structure
    ifgdir = os.path.join(procdir,'IFG')
    if not os.path.exists(ifgdir):
        os.mkdir(ifgdir)
    ifgthisdir = os.path.join(ifgdir,pair)
    if os.path.exists(ifgthisdir):
        shutil.rmtree(ifgthisdir)
    os.mkdir(ifgthisdir)
############################################################ remosaicking RSLCs (may not exist)
    for pomdate in [masterdate.strftime('%Y%m%d'),slavedate.strftime('%Y%m%d')]:
        if (not is_non_zero_file(os.path.join(procdir,'RSLC',pomdate,pomdate+'.rslc.par'))) \
            or (not os.path.exists(os.path.join(procdir,'RSLC',pomdate,pomdate+'.rslc'))):
            print('Regenerating mosaic for '+pomdate)
            rslctab=os.path.join(procdir,'tab',
                             pomdate+'R_tab')
            filename=os.path.join(procdir,'RSLC',pomdate,pomdate+'.rslc')
            logfile = os.path.join(procdir,'log',"mosaic_rslc_{0}.log".format(pomdate))
            #donno why but there is error here - glob.glob is not found!
            import glob
            iwfiles = glob.glob(os.path.join(procdir,'RSLC',pomdate,pomdate+'.IW*.rslc'))
            swathlist = []
            for iwfile in iwfiles:
                swathlist.append(iwfile.split('/')[-1].split('.')[1])
            swathlist.sort()
            rc, msg = make_SLC_tab(rslctab,filename,swathlist)
            if rc > 0:
                print('Something went wrong creating the slave resampled tab file...')
                return 1
            SLC_mosaic_S1_TOPS(rslctab,filename,rglks,azlks,logfile,mastertab)
############################################################ Create offsets
    offsetfile = os.path.join(ifgthisdir,pair+'.off')
    masterpar = os.path.join(procdir,'RSLC',
                         masterdate.strftime('%Y%m%d'),
                         masterdate.strftime('%Y%m%d')+'.rslc.par')
    slavepar = os.path.join(procdir,'RSLC',
                        slavedate.strftime('%Y%m%d'),
                        slavedate.strftime('%Y%m%d')+'.rslc.par')
    logfilename = os.path.join(procdir,'log','create_offset_{0}.log'.format(pair))
    
    if not create_offset(masterpar,slavepar,offsetfile,rglks,azlks,logfilename):
        print('\nERROR:', file=sys.stderr)
        print('\nSomething went wrong with creating the offset file {0}.'.format(offsetfile), file=sys.stderr)
        shutil.rmtree(ifgthisdir)
        return 1
############################################################ Est. Topographic phase
    origmasterpar = os.path.join(procdir,'RSLC',
                                 origmasterdate.strftime('%Y%m%d'),
                                 origmasterdate.strftime('%Y%m%d')+'.rslc.par')
    hgtfile = os.path.join(procdir,'geo',origmasterdate.strftime('%Y%m%d')+'.hgt')
    simfile = os.path.join(ifgthisdir,pair+'.sim_unw')
    logfilename = os.path.join(procdir,'log','phase_sim_orb_{0}.log'.format(pair))
    #if rglks != gc.rglks or azlks != gc.azlks:
    #    hgtfile = os.path.join(procdir,'geo',origmasterdate.strftime('%Y%m%d')+'.hgt.'+str(rglks)+'.'+str(azlks))
    #    if not os.path.exists(hgtfile):
    #        print('recomputing hgt file for custom looks')
    #        
    print('Estimating topographic phase...')
    if not phase_sim_orb(masterpar,slavepar,origmasterpar,offsetfile,hgtfile,simfile,logfilename):
        print('\nERROR:', file=sys.stderr)
        print('\nSomething went wrong with estimating the topographic phase estimation.', file=sys.stderr)
        shutil.rmtree(ifgthisdir)
        return 2

############################################################ Create interferogram
    print('Forming interferogram...')
    difffile = os.path.join(ifgthisdir,pair+'.diff')
    logfilename = os.path.join(procdir,'log','SLC_diff_intf_{0}.log'.format(pair))
    if not SLC_diff_intf(masterpar[:-4],slavepar[:-4],masterpar,slavepar,offsetfile,simfile,difffile,rglks,azlks,logfilename):
        print('\nERROR:', file=sys.stderr)
        print('\nSomething went wrong during the interferogram formation.', file=sys.stderr)
        print('try yourself by: SLC_diff_intf '+masterpar[:-4]+' '+slavepar[:-4]+' '+masterpar+' '+slavepar+' '+offsetfile+' '+simfile+' '+difffile+' '+str(rglks)+' '+str(azlks))
        shutil.rmtree(ifgthisdir)
        return 3

#create a ras file
    logfilename = os.path.join(procdir,'log','rasmph_pwr_{0}.log'.format(pair))
    if not rasmph_pwr(difffile,mastermli[:-4],str(width),logfilename):
        print('\nERROR:', file=sys.stderr)
        print('\nSomething went wrong during the rasmph step.', file=sys.stderr)
        print('\nTry yourself: rasmph_pwr '+difffile+' '+mastermli[:-4]+' '+str(width))
        shutil.rmtree(ifgthisdir)
        return 4

############################################################ Log IFG to Database
    if job_id != -1:
        print('About to insert new IFG product with job_id: %d' % job_id)
        lq.set_new_ifg_product(job_id, masterpar[:-4],slavepar[:-4], difffile)

############################################################ Estimate coherence
    mastermli = os.path.join(procdir,'RSLC',
                             masterdate.strftime('%Y%m%d'),
                             masterdate.strftime('%Y%m%d')+'.rslc.mli')
    slavemli = os.path.join(procdir,'RSLC',
                             slavedate.strftime('%Y%m%d'),
                             slavedate.strftime('%Y%m%d')+'.rslc.mli')
    #get baselines information to the qualityfile
    #templog = procdir+'/log/tmp_base.log'
    templog = procdir+'/log/tmp_base_'+masterdate.strftime('%Y%m%d')+'_'+slavedate.strftime('%Y%m%d')+'.log'
    pom = os.system('base_orbit '+mastermli+'.par '+slavemli+'.par - > '+templog)
    with open(qualityfile, "a") as myfile:
            myfile.write(grep1('perpendicular', templog))
            myfile.write(grep1('parallel', templog))
    #os.remove(procdir+'/log/tmp_base.log')
    os.remove(templog)
    
    cohfile = os.path.join(ifgdir,pair,pair+'.cc')
    logfilename  = os.path.join(procdir,'log','cc_wave_{0}.log'.format(pair))
    if not cc_wave(difffile,mastermli,slavemli,cohfile,width,logfilename):
        print('\nERROR:', file=sys.stderr)
        print('\nSomething went wrong during the coherence estimation.', file=sys.stderr)
        shutil.rmtree(ifgthisdir)
        return 5
    
    #create a coherence ras file
    logfilename = os.path.join(procdir,'log','rascc_{0}.log'.format(pair))
    if not rascc(cohfile,mastermli[:-4],str(width),1,logfilename):
        print('\nERROR:', file=sys.stderr)
        print('\nSomething went wrong during the coherence sunraster creation.', file=sys.stderr)
        shutil.rmtree(ifgthisdir)
        return 1


############################################################ Log coherence to database
    if job_id != -1:
        lq.set_new_coherence_product(job_id, masterpar[:-4], slavepar[:-4], difffile)

    return 0

################################################################################
#make baselines function
################################################################################
def make_baselines(procdir, masterdate, ziplistfile, lq):
    """ Use the gamma base_calc to generate a baseline list for interferogram selection/generation
    """
    # Filepaths:
    print('Estimating baeslines...')
    masterpar = os.path.join(procdir,'RSLC',
                             masterdate.strftime('%Y%m%d'),
                             masterdate.strftime('%Y%m%d')+'.rslc.par')
    logfilename = os.path.join(procdir,'log','base_calc_'+masterdate.strftime('%Y%m%d')+'.log')
    bperpfile = os.path.join(procdir,'bperp_'+masterdate.strftime('%Y%m%d')+'_file')
    itabfile = os.path.join(procdir,'itab_'+masterdate.strftime('%Y%m%d')+'')
    
    # Retrieve datelist (ie. based on ziplist)
    datelist = lq.get_zip_dates(lq.loadziplist(ziplistfile))
    datelist.sort()
    # Generate an SLC_tab
    slctab = os.path.join(procdir,'tab','base_calc_{0}_SLC_tab'.format(masterdate))
    try:
        with open(slctab, 'w') as f:
            for date in datelist:
                rslcf = os.path.join(procdir,'RSLC',
                                    date.strftime('%Y%m%d'),
                                    date.strftime('%Y%m%d')+'.rslc')
                if os.path.exists(rslcf) and os.path.exists(rslcf+'.par') and os.path.getsize(rslcf)>0 and os.path.getsize(rslcf+'.par') > 0:
                    f.write(rslcf+" "+rslcf+".par\n")
    except:
        print('\nERROR:', file=sys.stderr)
        print('\nFailed to write SLC_tab: '+slctab, file=sys.stderr)
        return 1
    
    if not base_calc(slctab,masterpar,bperpfile,itabfile,gc.Bp_min, gc.Bp_max, gc.Bt_min, gc.Bt_max, gc.nmax, logfilename):
        print('\nERROR:', file=sys.stderr)
        print('\nSomething went wrong with estimating baselines.', file=sys.stderr)
        return 1
    else:
        print('\nBaseline estimation complete. Written to: '+bperpfile, file=sys.stderr)

################################################################################
#make interferogram list
################################################################################
def make_interferograms_list(procdir, origmasterdate, bperpfile):
    print('Building interferograms based on list: '+bperpfile)
    
    if not( os.path.exists(bperpfile) and os.path.getsize(bperpfile)>0 ):
        print('\nERROR:', file=sys.stderr)
        print('\nCannot find baseline file!', file=sys.stderr)
        return 1
    else:
        with open(bperpfile, 'r') as f:
            for line in f:
                mdate = line.split()[1]
                sdate = line.split()[2]
                masterdate = dt.date(int(mdate[:4]),int(mdate[4:6]),int(mdate[6:8]))
                slavedate = dt.date(int(sdate[:4]),int(sdate[4:6]),int(sdate[6:8]))
                make_interferogram(origmasterdate,masterdate,slavedate,procdir, -1)

