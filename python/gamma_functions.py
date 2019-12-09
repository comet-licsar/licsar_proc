#!/usr/bin/env python
""" Wrapper and helper functions for the Gamma InSAR software package

=========
Changelog
=========
Sept 2016: Original implementation (KS)
Feb 2019: Update to use gamma/20181130 version (ML)
============
Contributors
============
Karsten Spaans, University of Leeds
Milan Lazecky, University of Leeds
"""

import os
import sys
import subprocess as subp
from glob import glob
import pdb
import global_config as gc


def rasrmg(uwfile,mlifile,width,reducfac,logfilename):
    """
    """
    rmgcall = ['rasrmg',uwfile,mlifile,str(width),'-','-','-',str(reducfac),str(reducfac)]
    with open(logfilename,'w') as f:
        try:
            rc = subp.check_call(rmgcall,stdout=f)
        except:
            print('Something went wrong during the unwrapped phase preview generation. Log file {0}'.format(logfilename))
            return False
    if rc != 0:
        print('Something went wrong during the unwrapped phase preview generation. Log file {0}'.format(logfilename))
        return False
    else:
        return True


def mcf(difffile,ccfile,mask,uwoutfile,width,logfilename):
    """
    """
    mcfcall = ['mcf',difffile,ccfile,mask,uwoutfile,str(width),'0','-','-','-','-','1','1']
    with open(logfilename,'w') as f:
        try:
            rc = subp.check_call(mcfcall,stdout=f)
        except:
            print('Something went wrong during the unwrapping. Log file {0}'.format(logfilename))
            return False

    if rc != 0:
        print('Something went wrong during the unwrapping. Log file {0}'.format(logfilename))
        return False
    else:
        return True

def look_vector(slcpar,offpar,dempar,dem,theta,phi,logfilename):
    """
    """
    lookcall = ['look_vector',slcpar,offpar,dempar,dem,theta,phi]
    with open(logfilename,'w') as f:
        try:
            rc = subp.check_call(lookcall,stdout=f)
        except:
            print('Something went wrong during the look vector calculation. Log file {0}'.format(logfilename))
            return False
    if rc != 0:
        print('Something went wrong during the look vector calculation. Log file {0}'.format(logfilename))
        return False
    else:
        return True

def adf(difffile,filtout,ccout,width,alpha,winsize,logfilename):
    """
    """
    adfcall = ['adf',difffile,filtout,ccout,str(width),str(alpha),str(winsize)]
    with open(logfilename,'w') as f:
        try:
            rc = subp.check_call(adfcall,stdout=f)
        except:
            print('Something went wrong during the Goldstein filtering. Log file {0}'.format(logfilename))
            return False

    if rc != 0:
        print('Something went wrong during the Goldstein filtering. Log file {0}'.format(logfilename))
        return False
    else:
        return True


def rascc(cohfile,mastermli,width,reducfac,logfilename):
    """
    """
    rascall = ['rascc',cohfile,mastermli,str(width),'-','-','-',str(reducfac),str(reducfac),'0.0','1.0','0.8','0.35']
    with open(logfilename,'w') as f:
        try:
            rc = subp.check_call(rascall,stdout=f)
        except:
            print('Something went wrong during the coherence sunraster creation. Log file {0}'.format(logfilename))
            return False

    if rc != 0:
        print('Something went wrong during the coherence sunraster creation. Log file {0}'.format(logfilename))
        return False
    else:
        return True

def cc_wave(difffile,mastermli,slavemli,cohfile,width,logfilename):
    """
        Makes coherence
    """
    cccall = ['cc_wave',difffile,mastermli,slavemli,cohfile,str(width)]
    with open(logfilename,'w') as f:
        try:
            rc = subp.check_call(cccall,stdout=f)
        except:
            print('Something went wrong during the coherence estimation. Log file {0}'.format(logfilename))
            return False

    if rc != 0:
        print('Something went wrong during the coherence estimation. Log file {0}'.format(logfilename))
        return False
    else:
        return True

def SLC_diff_intf(masterslc,slaveslc,masterpar,slavepar,offfile,simfile,difffile,rglks,azlks,logfilename):
    """
        Makes interferogram
    """
    intcall = ['SLC_diff_intf',masterslc,slaveslc,masterpar,slavepar,offfile,simfile,difffile,str(rglks),str(azlks),'0','0','0.2','1','1']

    with open(logfilename,'w') as f:
        try:
            rc = subp.check_call(intcall,stdout=f)
        except:
            print('Something went wrong during the interferogram formation. Log file {0}'.format(logfilename))
            return False
    if rc != 0:
        print('Something went wrong during the interferogram formation. Log file {0}'.format(logfilename))
        return False
    else:
        return True


def phase_sim_orb(masterpar,slavepar,origmasterpar,offfile,hgtfile,simfile,logfilename):
    """
    """
    simcall = ['phase_sim_orb',masterpar,slavepar,offfile,hgtfile,
               simfile,origmasterpar,'-','-','1','1']
    with open(logfilename,'w') as f:
        try:
            rc = subp.check_call(simcall,stdout=f)
        except:
            print('Something went wrong during the simulating the topographic phase. Log file {0}'.format(logfilename))
            return False
    if rc != 0:
        print('Something went wrong during the simulating the topographic phase. Log file {0}'.format(logfilename))
        return False
    else:
        return True

def S1_coreg_overlap(mastertab,slavetab,pair,offfile,offrefine,auxtab,logfilename):
    """ Estimate azimuth offset between master and slave image used ESD
    """
    if auxtab == '':
        overlapcall = ['S1_coreg_overlap',mastertab,slavetab,pair,offfile,offrefine,
                   '0.8','0.01','0.8','1']
    else:
        overlapcall = ['S1_coreg_overlap',mastertab,slavetab,pair,offfile,offrefine,
                   '0.8','0.01','0.8','1',auxtab]
    with open(logfilename,'w') as f:
        try:
            rc = subp.check_call(overlapcall,stdout=f)
        except:
            print('Something went wrong during the spectral diversity offset estimation. Log file {0}'.format(logfilename))
            return False

    if rc != 0:
        print('Something went wrong during the spectral diversity offset estimation. Log file {0}'.format(logfilename))
        return False
    else:
        return True


def S1_coreg_subswath_overlap(mastertab,slavetab,pair,offfile,offrefine,auxtab,logfilename):
    """ Estimate coregistration offset based on the subswath overlap
    """
    if auxtab == '':
        overlapcall = ['S1_coreg_subswath_overlap',mastertab,slavetab,pair,offfile,offrefine,
                   '0.8','0.01','0.8','1']
    else:
        overlapcall = ['S1_coreg_subswath_overlap',mastertab,slavetab,pair,offfile,offrefine,
                   '0.8','0.01','0.8','1',auxtab]
    with open(logfilename,'w') as f:
        try:
            rc = subp.check_call(overlapcall,stdout=f)
        except:
            print('Something went wrong during the subswath overlap function. Log file {0}'.format(logfilename))
            return False

    if rc != 0:
        print('Something went wrong during the subswath overlap function. Log file {0}'.format(logfilename))
        return False
    else:
        return True


def create_offset(masterpar,slavepar,offfile,rglks,azlks,logfilename):
    """
    """
    offsetcall = ['create_offset',masterpar,slavepar,offfile,'1',str(rglks),str(azlks),'0']
    with open(logfilename,'w') as f:
        try:
            rc = subp.check_call(offsetcall,stdout=f)
        except:
            print('Something went wrong creating the offset file. Log file {0}'.format(logfilename))
            return False
    if rc != 0:
        print('Something went wrong creating the offset file. Log file {0}'.format(logfilename))
        return False
    return True


def SLC_interp_lt_S1_TOPS(slavetab,slavepar,mastertab,masterpar,lut,mastermlipar,slavemlipar,offpar,slaveRtab,slaveRslc,slaveRpar,logfilename):
    """ Resampling of SLC based on lookup table and offset file
    """
    interpcall = ['SLC_interp_lt_ScanSAR',slavetab,slavepar,mastertab,
                  masterpar,lut,mastermlipar,slavemlipar,offpar,
                  slaveRtab,slaveRslc,slaveRpar]
    with open(logfilename,'w') as f:
        try:
            rc = subp.check_call(interpcall,stdout=f)
        except:
            print('Something went wrong resampling image. Log file {0}'.format(logfilename))
            return False
    if rc != 0:
        print('Something went wrong resampling image. Log file {0}'.format(logfilename))
        return False
    return True

def SLC_interp(slaveslc,masterpar,slavepar, offpar, slaveRslc, slaveRpar,logfilename):
    """ Resampling of SLC based on lookup table and offset file
    """
    interpcall = ['SLC_interp',slaveslc, masterpar,slavepar,
                  offpar, slaveRslc, slaveRpar]
    with open(logfilename,'w') as f:
        try:
            rc = subp.check_call(interpcall,stdout=f)
        except:
            print('Something went wrong resampling image. Log file {0}'.format(logfilename))
            return False
    if rc != 0:
        print('Something went wrong resampling image. Log file {0}'.format(logfilename))
        return False
    return True

def rdc_trans(mastermlipar,demhgt,slavemlipar,lut,logfilename):
    """
    """
    rdccall = ['rdc_trans',mastermlipar,demhgt,slavemlipar,lut]
    with open(logfilename,'w') as f:
        try:
            rc = subp.check_call(rdccall,stdout=f)
        except:
            print('Something went wrong deriving the lookup table from master to slave geometry. Log file {0}'.format(logfilename))
            return False
    if rc != 0:
        print('Something went wrong deriving the lookup table from master to slave geometry. Log file {0}'.format(logfilename))
        return False
    else:
        return True

def rasmph_pwr(difffile,mastermli,width,logfilename):
    """
    """
    rascall = ['rasmph_pwr',difffile,mastermli,str(width),'-','-','-','1','1','0.8','0.35']
    with open(logfilename,'w') as f:
        try:
            rc = subp.check_call(rascall,stdout=f)
        except:
            print('Something went wrong creating sunraster file. Log file {0}'.format(logfilename))
            return False
    if rc != 0:
        print('Something went wrong creating sunraster file. Log file {0}'.format(logfilename))
        return False
    else:
        return True


def rashgt(hgt,pwr,width,starthgt,startpwr,nlines,pixavr,pixavaz,mcycle,logfilename):
    """
    """
    rascall = ['rashgt',hgt,pwr,str(width),starthgt,startpwr,str(nlines),pixavr,pixavaz,mcycle]
    with open(logfilename,'w') as f:
        try:
            rc = subp.check_call(rascall,stdout=f)
        except:
            print('Something went wrong creating sunraster file. Log file {0}'.format(logfilename))
            return False
    if rc != 0:
        print('Something went wrong creating sunraster file. Log file {0}'.format(logfilename))
        return False
    else:
        return True

def geocode(gcmap,datain,widthin,dataout,widthout,nlinesout,interpmode,dtype,logfilename):
    """
    """
    geocodecall = ['geocode',gcmap,datain,str(widthin),dataout,str(widthout),
                   str(nlinesout),interpmode,dtype]
    with open(logfilename,'w') as f:
        try:
            rc = subp.check_call(geocodecall,stdout=f)
        except:
            print('Something went wrong creating master height file. Log file {0}'.format(logfilename))
            return False
    if rc != 0:
        print('Something went wrong creating master height file. Log file {0}'.format(logfilename))
        return False
    else:
        return True
def S1_OPOD_vec(parFile,orbitFile,logfilename):
    """
    Correct the slc parameter files
    """
    OPOD_call = ['S1_OPOD_vec',parFile,orbitFile]
    with open(logfilename,'w') as f:
        try:
            rc = subp.check_call(OPOD_call,stdout=f)
        except:
            print('Something went wrong corecting {0} with orbit {1} - see log {2}'.format(parFile,orbitFile,logfilename))
            return False
      # S1_OPOD_vec ${img} $orbit_file 
    return True 

def geocode_back(datain,widthin,gcmap,dataout,widthout,nlinesout,interpmode,dtype,logfilename):
    """
    """
    backcall = ['geocode_back',datain,str(widthin),gcmap,dataout,str(widthout),
                str(nlinesout),str(interpmode),dtype]
    with open(logfilename,'w') as f:
        try:
            rc = subp.check_call(backcall,stdout=f)
        except:
            print('Something went wrong creating the DEM mli. Log file {0}'.format(logfilename))
            return False
    if rc != 0:
        print('Something went wrong creating the DEM mli. Log file {0}'.format(logfilename))
        return False
    else:
        return True

def offset_fitm(offs,ccp,diffpar,coffs,coffsets,thres,npoly,logfilename):
    """
    """
    fitcall = ['offset_fitm',offs,ccp,diffpar,coffs,coffsets,str(thres),str(npoly)]
    with open(logfilename,'w') as f:
        try:
            rc = subp.check_call(fitcall,stdout=f)
        except:
            print('Something went wrong during the offset function fitting. Log file {0}'.format(logfilename))
            return False
    if rc != 0:
        print('Something went wrong during the offset function fitting. Log file {0}'.format(logfilename))
        return False
    else:
        return True

def offset_fit(p,doffset,logfilename,acqMode='iw'):
    """
    """
    if acqMode == 'iw':
        thres = 0.2
        npoly = 1
    elif acqMode == 'sm':
        thres = 0.15
        npoly = 3

    fitcall = ['offset_fit',p+'.offs',p+'.snr',doffset,'-','-',str(thres),str(npoly),'0']
    with open(logfilename,'w') as f:
        try:
            rc = subp.check_call(fitcall,stdout=f)
        except:
            print('Something went wrong during the offset function fitting. Log file {0}'.format(logfilename))
            return False
    if rc != 0:
        print('Something went wrong during the offset function fitting. Log file {0}'.format(logfilename))
        return False
    else:
        return True

def set_value(parin,parout,keyw,newval):
    """
    """
    setcall = ['set_value',parin,parout,keyw,str(newval)]
    rc = subp.check_call(setcall)
    if rc != 0:
        print('Something went wrong during setting value.')
        return False
    else:
        return True

def offset_pwrm(mli1,mli2,diffpar,offfile,ccpfile,rwin,azwin,offsets,n_ovr,nr,naz,thres,logfilename):
    """
    """
    pwrmcall = ['offset_pwrm',mli1,mli2,diffpar,offfile,ccpfile,str(rwin),str(azwin),offsets,str(n_ovr),str(nr),str(naz),str(thres)]
    with open(logfilename,'w') as f:
        try:
            rc = subp.check_call(pwrmcall,stdout=f)
        except:
            print('Something went wrong during the cross correlation offset estimation. Log file {0}'.format(logfilename))
            return False
    if rc != 0:
        print('Something went wrong during the cross correlation offset estimation. Log file {0}'.format(logfilename))
        return False
    else:
        return True

def offset_pwr_tracking(slc1,rslc,slc1_par,rslc_par,dofffile,p,rstep,azstep,logfilename):
    """
    """
    pwrcall = ['offset_pwr_tracking',slc1,rslc,slc1_par,rslc_par,dofffile,p+'.offs',p+'.snr','128','64','-','1','0.2',str(rstep),str(azstep),'-','-','-','-','-','-','1']
    with open(logfilename,'w') as f:
        try:
            rc = subp.check_call(pwrcall,stdout=f)
        except:
            print('Something went wrong during the cross correlation offset tracking. Log file {0}'.format(logfilename))
            return False
    if rc != 0:
        print('Something went wrong during the cross correlation offset tracking. Log file {0}'.format(logfilename))
        return False
    else:
        return True

def offset_pwr(slc1,slc2,slc1_par,slc2_par,offfile,p,logfilename):
    """
    """
    pwrcall = ['offset_pwr',slc1,slc2,slc1_par,slc2_par,offfile,p+'.offs',p+'.snr','32','32','-','1','-','-','0.15']
    with open(logfilename,'w') as f:
        try:
            rc = subp.check_call(pwrcall,stdout=f)
        except:
            print('Something went wrong during the cross correlation offset estimation. Log file {0}'.format(logfilename))
            return False
    if rc != 0:
        print('Something went wrong during the cross correlation offset estimation. Log file {0}'.format(logfilename))
        return False
    else:
        return True

def create_diff_par(par1,par2,diffpar,partype,iflag,logfilename):
    """
    """
    parcall = ['create_diff_par',par1,par2,diffpar,partype,iflag]
    with open(logfilename,'w') as f:
        try:
            rc = subp.check_call(parcall,stdout=f)
        except:
            print('Something went wrong during the diff parameter file creation. Log file {0}'.format(logfilename))
            return False
    if rc != 0:
        print('Something went wrong during the diff parameter file creation. Log file {0}'.format(logfilename))
        return False
    else:
        return True


def pixel_area(mlipar,dem,lut,lsmap,incmap,pixsigma,pixgamma,logfilename):
    """
    """
    pixelcall = ['pixel_area',mlipar,dem+'_par',dem,lut,lsmap,incmap,pixsigma,pixgamma]
    with open(logfilename,'w') as f:
        try:
            rc = subp.check_call(pixelcall,stdout=f)
        except:
            print('Something went wrong during the DEM amplitude simulation. Log file {0}'.format(logfilename))
            return False
    if rc != 0:
        print('Something went wrong during the DEM amplitude simulation. Log file {0}'.format(logfilename))
        return False
    else:
        return True

def poly_math(mastermli,mastermliaz,mastermliwidth,masterazpolyout,logfilename):
    """
    """
    polycall = ['poly_math',mastermli,str(mastermliaz),str(mastermliwidth),masterazpolyout,'-','1','0.0','1.0']
    with open(logfilename,'w') as f:
        try:
            rc = subp.check_call(polycall,stdout=f)
        except:
            print('Something went wrong during the poly_math. Log file {0}'.format(logfilename))
            return False
    if rc != 0:
        print('Something went wrong during the poly_math. Log file {0}'.format(logfilename))
        return False
    else:
        return True



def mask_class(mazovr,lut,lutazovr,logfilename):
    """
    """
    mask_call=['mask_class',mazovr,lut,lutazovr,'1','1','1','1','0','0.0','0.0']
    with open(logfilename,'w') as f:
        try:
            rc = subp.check_call(mask_call,stdout=f)
        except:
            print('Something went wrong during the mask_class. Log file {0}'.format(logfilename))
            return False
    if rc != 0:
        print('Something went wrong during the mask_class. Log file {0}'.format(logfilename))
        return False
    else:
        return True


def raspwr(image,wid,pixavr,pixavaz,rasout,logfilename):
    """
    """
    rascall = ['raspwr',image,str(wid),'1','0',str(pixavr),str(pixavaz),'1.0','0.35','1',rasout]
    with open(logfilename,'w') as f:
        try:
            rc = subp.check_call(rascall,stdout=f)
        except:
            print('Something went wrong during the raspwr. Log file {0}'.format(logfilename))
            return False
    if rc != 0:
        print('Something went wrong during the raspwr. Log file {0}'.format(logfilename))
        return False
    else:
        return True

def S1_poly_overlap(mastertab,rlk,azlk,polyout,overlap_type,logfilename):
    """
    """
    polycall = ['S1_poly_overlap',mastertab,str(rlk),str(azlk),polyout,str(overlap_type)]
    with open(logfilename,'w') as f:
        try:
            rc = subp.check_call(polycall,stdout=f)
        except:
            print('Something went wrong during the S1_poly_overlap. Log file {0}'.format(logfilename))
            return False
    if rc != 0:
        print('Something went wrong during the S1_poly_overlap. Log file {0}'.format(logfilename))
        return False
    else:
        return True

def gc_map_fine(gcin,width,diffpar,gcout,refflag,logfilename):
    """
    """
    mapcall = ['gc_map_fine',gcin,str(width),diffpar,gcout,'1']
    with open(logfilename,'w') as f:
        try:
            rc = subp.check_call(mapcall,stdout=f)
        except:
            print('Something went refining the lookup table. Log file {0}'.format(logfilename))
            return False
    if rc != 0:
        print('Something went refining the lookup table. Log file {0}'.format(logfilename))
        return False
    else:
        return True

def gc_map(mlipar,offpar,dem,demseg,lut,latovr,lonovr,simsar,u,v,inc,psi,pix,lsmap,frame,lsmode,r_ovr,logfilename):
    """
    """
    mapcall = ['gc_map',mlipar,offpar,dem+'_par',dem,demseg+'_par',demseg,lut,
               latovr,lonovr,simsar,u,v,inc,psi,pix,lsmap,frame,lsmode,str(r_ovr)]
    with open(logfilename,'w') as f:
        try:
            rc = subp.check_call(mapcall,stdout=f)
        except:
            print('Something went wrong geocoding the master mli image. Log file {0}'.format(logfilename))
            return False
    if rc != 0:
        print('Something went wrong geocoding the master mli image. Log file {0}'.format(logfilename))
        return False
    else:
        return True

def SLC_mosaic_S1_TOPS(tabname,outfilename,rglks,azlks,logfilename,mastertab=''):
    """
    Include mastertab when mosaicing a resampled SLC
    """
    if mastertab == '':
        mosaiccall = ['SLC_mosaic_S1_TOPS',tabname,outfilename,
                      outfilename+'.par',str(rglks),str(azlks)]
    else:
        mosaiccall = ['SLC_mosaic_S1_TOPS',tabname,outfilename,
                      outfilename+'.par',str(rglks),str(azlks),
                      '0',mastertab]
    with open(logfilename,'w') as f:
        try:
            rc = subp.check_call(mosaiccall,stdout=f)
        except:
            print('Something went wrong mosaicing file {0}. Log file {1}'.format(outfilename,logfilename))
            return False
    if rc != 0:
        print('Something went wrong mosaicing file {0}. Log file {1}'.format(outfilename, logfilename))
        return False
    else:
        return True

def mosaic_TOPS(procdir,imdir,imdate,swathlist):
    """ Mosaic the subswaths together

    In:
        procdir     path to main processing directory
        imdir       path to procdir/SLC/imdate directory
        imdate      datetime.date object with image date
        swathlist   list of swaths to be mosaiced
    Out:
        boolean True if successful, False if not

    """
    # Check if more than 1 subswath in swathlist
    if len(swathlist) < 2 and 'IW' in swathlist[0]:
        return True
    print('Mosaicing subswaths...')
    tabname = os.path.join(procdir,'tab',imdate.strftime('%Y%m%d')+'_tab')
    filename = os.path.join(imdir,imdate.strftime('%Y%m%d'))
    make_SLC_tab(tabname,filename,swathlist)
    logfilename = os.path.join(procdir,'log','mosaic_TOPS_{0}.log'.format(imdate.strftime('%Y%m%d')))
    mosaiccall = ['SLC_mosaic_S1_TOPS',tabname,os.path.join(imdir,imdate.strftime('%Y%m%d')+'.slc'),
                  os.path.join(imdir,imdate.strftime('%Y%m%d')+'.slc.par'),str(rglks),str(azlks)]
    with open(logfilename,'w') as f:
        try:
            rc = subp.check_call(mosaiccall,stdout=f)
        except:
            print('Something went wrong mosaicing swaths {0}. Log file {1}'.format(' '.join([sw for sw in swathlist]),logfilename))
            #return False
    #Perhaps something extra in stdout after gamma update? it seems working but rc is not 0
    if rc != 0:
        print('Something went wrong mosaicing swaths {0}. Log file {1}'.format(' '.join([sw for sw in swathlist]),logfilename))
        #return False

    # just trying to continue - seems that mosaics exist even though an error is reported
    logfilename = os.path.join(procdir,'log','multilookSLC_{0}.log'.format(imdate.strftime('%Y%m%d')))
    multicall = 'multilookSLC {0} {1} {2} 1 {3} &> {4}'.format(imdate.strftime('%Y%m%d'),rglks,azlks,imdir,logfilename)
    rc = os.system(multicall)
    if rc == 0:
        return True
    else:
        print('Something went wrong multilooking the SLC. Continuing with next image.')
        return False

def SLC_copy_S1_TOPS(SLColdtab,SLCnewtab,bursttab,imdir,procdir,logfile):
    """Extract bursts from subswath file and copies to new file

    In:
        SLColdtab    SLC tab of original file containing bursts
        SLCnewtab    SLC tab of destination file
    bursttab     burst tab
        imdir        path to SLC imdate directory *** CAN BE REMOVED?? ***
        procdir      path to main processing directory
        logfile      path and filename of logfile where output is written to
    Out:
        0 if programme finished successfully, 1 if not
    """
    copycall = ['SLC_copy_ScanSAR',SLColdtab,SLCnewtab,bursttab]
    try:
        with open(logfile,'w') as lf:
            rc = subp.check_call(copycall,stdout=lf)
    except:
        print('Could not copy bursts from file in SLC tab file {0}. Log file {1}.'.format(SLCnewtab,logfile))
        return 1
    if rc != 0:
        print('Could not copy bursts from file in SLC tab file {0}. Log file {1\}.'.format(SLCnewtab,logfile))

    return 0

def par_S1_SLC(tiff,annot,calib,noise,slc,logfile,acqMode='iw'):
    """
    """
    if acqMode == 'iw':
        topsfile = slc+'.TOPS_par'
        dtype = 1 
    elif acqMode == 'sm':
        topsfile = '-'
        dtype = 1
    parcall = ['par_S1_SLC',tiff,annot,calib,noise,
                       slc+'.par',slc,topsfile,str(dtype),'60']
    try:
        with open(logfile,'w') as lf:
            rc = subp.check_call(parcall,stdout=lf,stderr=lf)
    except:
        print('There was a problem converting file {0} into gamma format. Log file {1}. Skipping date.'.format(tiff,logfile))
        return False
    if rc != 0:
        print('There was a problem converting file {0} into gamma format. Log file {1}. Skipping date.'.format(tiff,logfile))
        return False
    return True

def make_SLC_tab(tabname,filename,swathlist):
    """Makes SLC_tab input file for Gamma

    In:
        tabname     Name of SLC_tab file
        filename    Name of image file
        swathlist   List of subswaths to be included in tab
    Out:
        0 if successful, 1 and error string if not
    """
    filebasename = '.'.join(filename.split('.')[:-1])
    fileextension = filename.split('.')[-1]
    try:
        with open(tabname,'w') as f:
            for s in swathlist:
                f.write('{0}.{1}.{2} {0}.{1}.{2}.par {0}.{1}.{2}.TOPS_par\n'.format(filebasename,s,fileextension))
    except:
        e = sys.exc_info()[0]
        return 1,e
    else:
        return 0,''

def make_burst_tab(tabname,startbursts,endbursts):
    """Makes burst_tab input file for Gamma

    In:
        tabname     Name of burst_tab file
    startbursts  Start burst for each swath to extract, can be scalar or list
    endbursts    end burst for each swath to extract, can be scalar or list
    Out:
        0 if successful, 1 and error string if not
    """
    try:
        with open(tabname,'w') as f:
            if isinstance(startbursts,list):
                for sb,eb in map(None,startbursts,endbursts):
                    f.write('{0} {1}\n'.format(sb,eb))
            else:
                f.write('{0} {1}\n'.format(startbursts,endbursts))
    except:
        e = sys.exc_info()[0]
        return 1,e
    else:
        return 0,''

def make_mosaic_tab(tabname,filenames,swathlist):
    """Makes SLCmosaic_tab input file for Gamma

    In:
        tabname     Name of SLC_tab file
        filenames   Name of image files
    Out:
        0 if successful, 1 and error string if not
    """
    try:
        with open(tabname,'w') as f:
            for filename in filenames:
                filebasename = '.'.join(filename.split('.')[:-1])
                fileextension = filename.split('.')[-1]
                f.write('{0}.{1} {0}.{1}.par \n'.format(filebasename,fileextension))
    except:
        e = sys.exc_info()[0]
        return 1,e
    else:
        return 0,''

def SLC_cat_S1_TOPS(tab1,tab2,tab3,logfile):
    """Concatenates two SLC files together

    In:
        tab1       SLC tab of first file to be concatenated
        tab2       SLC tab of second file to be concatenated
        tab3       SLC tab of destination file
    Out:
        Boolean True if successful, False if not
    """

    catcall = ['SLC_cat_ScanSAR',tab1,tab2,tab3]
    with open(logfile,'w') as lf:
        try:
            rc = subp.check_call(catcall,stdout=lf)
        except:
            print('Could not concatenate files in tabs {0} and {1}. Log file {2}. Continuing with next acquisition date.'.format(tab1,tab2,logfile))
            return False
    if rc != 0:
        print('Could not concatenate files in tabs {0} and {1}. Log file {2}. Continuing with next acquisition date.'.format(tab1,tab2,logfile))
        return False
    return True

def SLC_phase_shift(SLC, SLC_par, SLC2, SLC2_par, ph_shift, logfile):
    """Correct phase shift
    
    Params:
        SLC       SLC data file to shift
        SLC_par   SLC parameter file
        SLC2      output SLC
        SLC2_par  output SLC parameter file
        ph_shift  phase shift to add to SLC phase (radians)

    NOTE: Used to apply a constant phase shift of -1.25 radians to Sentinel-1 TOPS SLC data
        from swath IW1 acquired up to 10-Mar-2015.
    """
    shiftcall = ['SLC_phase_shift',SLC,SLC_par,SLC2,SLC2_par,str(ph_shift)]
    with open(logfile,'w') as lf:
        try:
            rc = subp.check_call(shiftcall,stdout=lf)
        except:
            print('Could not shift phase of SLC file {0}. Log file {1}.'.format(SLC,logfile))
            return False
    if rc != 0:
        print('Could not shift phase of SLC file {0}. Log file {1}.'.format(SLC,logfile))
        return False
    return True
    
def base_calc(SLC_tab, SLC_par, bperp_file, itab, bperp_min, bperp_max, delta_T_min, delta_T_max, delta_n_max, logfile):
    """Creates baseline list for the RSLCs in SLC_tab
    In:
         SLC_tab      SLC tab
         SLC_par      SLC par file
         bperp_file   output bperp file
         itab
    Out:
        Boolean True if successful, False if not
    """

    catcall = ['base_calc', SLC_tab, SLC_par, bperp_file, itab, '1', '1', str(bperp_min), str(bperp_max), str(delta_T_min), str(delta_T_max), str(delta_n_max)]
    with open(logfile,'w') as lf:
        try:
            rc = subp.check_call(catcall,stdout=lf)
        except:
            print('Could make baseline list from files in tab: {0} and slc.par: {1}. Log file {2}. '.format(SLC_tab,SLC_par,logfile))
            print('Tried command:')
            print(' '.join(catcall))
            return False
    if rc != 0:
        print('Could make baseline list from files in tab: {0} and slc.par: {1}. Log file {2}. '.format(SLC_tab,SLC_par,logfile))
        print('Tried command:')
        print(' '.join(catcall))
        return False
    return True

def init_offset_orbit(masterpar,slavepar,offfile,logfilename):
    """
    """
    initofforbitcall = ['init_offset_orbit',masterpar,slavepar,offfile]
    with open(logfilename,'w') as f:
        try:
            rc = subp.check_call(initofforbitcall,stdout=f)
        except:
            print('Something went wrong estimating offset using orbit information. Log file {0}'.format(logfilename))
            return False
    if rc != 0:
        print('Something went wrong estimating offset using orbit information. Log file {0}'.format(logfilename))
        return False
    return True

def init_offset(masterfile,slavefile,masterpar,slavepar,offfile,rglks,azlks,logfilename):
    """
    """
    initoffcall = ['init_offset',masterfile,slavefile,masterpar,slavepar,offfile,str(rglks),str(azlks)]
    with open(logfilename,'w') as f:
        try:
            rc = subp.check_call(initoffcall,stdout=f)
        except:
            print('Something went wrong estimating offset using orbit information. Log file {0}'.format(logfilename))
            return False
    if rc != 0:
        print('Something went wrong estimating offset using orbit information. Log file {0}'.format(logfilename))
        return False
    return True

def multi_look(slc, slcpar, mli, mlipar, rlks, azlks, logfilename):
    """
    """
    mlcall = ['multi_look', slc, slcpar, mli, mlipar, str(rlks), str(azlks),'-','-','1e-6']
    with open(logfilename,'w') as f:
        try:
            rc = subp.check_call(mlcall,stdout=f)
        except:
            print('Something went wrong during multilooking. Log file {0}'.format(logfilename))
            return False
    if rc != 0:
        print('Something went wrong during multilooking. Log file {0}'.format(logfilename))
        return False
    return True



