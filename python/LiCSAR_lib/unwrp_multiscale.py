################################################################################
#Imports
################################################################################
import numpy as np
import os, glob
import gdal
import subprocess
import xarray as xr
import rioxarray

import cv2
from scipy import interpolate
from scipy import ndimage
import time
import matplotlib.pyplot as plt
xr.set_options(keep_attrs=True)
from LiCSAR_lib.LiCSAR_misc import *

from scipy.ndimage import gaussian_filter
from scipy.ndimage import generic_filter
from scipy import stats

from astropy.convolution import Gaussian2DKernel, interpolate_replace_nans, convolve_fft
from sklearn.linear_model import HuberRegressor

import shutil

#set dask client to use only one worker...
#from dask.distributed import Client
#client = Client(n_workers=1)

################################################################################

def interpolate_nans(array, method='cubic'):
    array = np.ma.masked_invalid(array)
    x = np.arange(0, array.shape[1])
    y = np.arange(0, array.shape[0])
    xx, yy = np.meshgrid(x, y)
    x1 = xx[~array.mask]
    y1 = yy[~array.mask]
    newarr = array[~array.mask]
    GD1 = interpolate.griddata((x1, y1), newarr.ravel(),(xx, yy),method=method) #, fill_value=0)
    GD1 = np.array(GD1)
    return GD1


def runcmd(cmd, printcmd = True):
    if printcmd:
        print(cmd)
        rc = os.system(cmd)
    else:
        #with nostdout():
        rc = os.system(cmd+' >/dev/null 2>/dev/null')
    #rc = os.system(cmd, shell=True, text=True)
    if not rc == 0:
        print('WARNING - command did not exit as OK')


def get_width(infile):
    #should work either with grd, nc or geotiffs
    batcmd="gmt grdinfo {} | grep n_columns".format(infile)
    result = subprocess.check_output(batcmd, shell=True, text=True)
    width = result.split()[-1]
    return width


def convert_tif2nc(tif, nc, coh=False):
    if coh:
        runcmd('gmt grdmath {0} 255 DIV = {1}'.format(tif, nc)) #out in px reg
    else:
        runcmd('gmt grdconvert {0} {1}; gmt grdedit {1} -T'.format(tif, nc))


def convert_nc2tif(nc, tif): #, coh=False):
    #if coh:
    #    runcmd('gmt grdmath {0} 255 DIV = {1}'.format(tif, nc)) #out in px reg
    #else:
    runcmd('gmt grdconvert {0} {1}; gmt grdedit {1} -T'.format(nc, tif))


def magpha2RI_array(mag, pha):
    R = np.cos(pha) * mag
    I = np.sin(pha) * mag
    out = R + 1j*I
    return out


def coh_from_cpx(cpxvalues, winsize = 3):
    #outcoh = cpx.copy(deep=True)
    realvar = ndimage.generic_filter(np.real(cpxvalues), np.var, size=winsize)
    imagvar = ndimage.generic_filter(np.imag(cpxvalues), np.var, size=winsize)
    fullvar = realvar + imagvar
    outcoh = 1/np.sqrt(1+2*winsize*winsize*fullvar)
    return outcoh


def filter_cpx_gauss(ifg_ml, sigma = 2, trunc = 4):
    R = np.real(ifg_ml.cpx.values)
    I = np.imag(ifg_ml.cpx.values)
    #
    gR = filter_nan_gaussian_conserving(R, sigma=sigma, trunc=trunc)
    gI = filter_nan_gaussian_conserving(I, sigma=sigma, trunc=trunc)
    #
    gauss_cpx = gR + 1j*gI
    gauss_xr = ifg_ml['cpx'].copy(deep=True)
    #
    gauss_xr.values = gauss_cpx
    #gauss_xr['pha'].values = np.angle(gauss_cpx)
    return gauss_xr



def filter_nan_gaussian_conserving(arr, sigma, trunc):
    """Apply a gaussian filter to an array with nans.
 thanks to:
 https://stackoverflow.com/questions/18697532/gaussian-filtering-a-image-with-nan-in-python
    Intensity is only shifted between not-nan pixels and is hence conserved.
    The intensity redistribution with respect to each single point
    is done by the weights of available pixels according
    to a gaussian distribution.
    All nans in arr, stay nans in gauss.
    """
    nan_msk = np.isnan(arr)
    loss = np.zeros(arr.shape)
    loss[nan_msk] = 1
    loss = ndimage.gaussian_filter(
            loss, sigma=sigma, truncate=trunc, mode='constant', cval=1)
    gauss = arr.copy()
    gauss[nan_msk] = 0
    gauss = ndimage.gaussian_filter(
            gauss, sigma=sigma, truncate=trunc, mode='constant', cval=0)
    gauss[nan_msk] = np.nan
    gauss += loss * arr
    return gauss


def modal(P):
    """We receive P[0]..P[8] with the pixels in the 3x3 surrounding window"""
    #if np.nanmin(P) == 0:
    out = stats.mode(P).mode[0]
    #else:
    #    out = 1
    return out


def filter_mask_modal(ifg_ml, masklayer = 'cohmask', outmask = 'cohmask2', winlen = 8):
    amask = ifg_ml[masklayer].where(ifg_ml.mask > 0).values
    #two iterations of modal filter - win 5 seemed optimal
    amask2 = generic_filter(amask, modal, (winlen, winlen))
    amask2 = amask2 * ifg_ml[masklayer].values
    #winlen = 5
    amask2 = generic_filter(amask2, modal, (winlen, winlen))
    ifg_ml[outmask] = ifg_ml[masklayer]
    ifg_ml[outmask].values = amask2 * ifg_ml[masklayer].values
    ifg_ml[outmask] = ifg_ml[outmask].fillna(0)
    return ifg_ml


def filter_modal(cohmask, winlen = 8):
    amask = ifg_ml.cohmask.where(ifg_ml.mask > 0).values
    #two iterations of modal filter - win 5 seemed optimal
    amask2 = generic_filter(amask, modal, (winlen, winlen))
    amask2 = amask2 * ifg_ml.cohmask.values
    #winlen = 5
    amask2 = generic_filter(amask2, modal, (winlen, winlen))
    ifg_ml['cohmask2'] = ifg_ml['cohmask']
    ifg_ml['cohmask2'].values = amask2 * ifg_ml.cohmask.values
    return ifg_ml


def create_preview_bin(binfile, width, ftype = 'unw'):
    # possible ftypes: unw, pha, coh, mag
    r=''
    if ftype == 'unw':
        q='normal'
        f='r4'
        c='jet'
        r='-r -20,20'
    elif ftype == 'pha':
        q='phase'
        f='cr4'
        c='jet'
    elif ftype == 'coh':
        q='normal'
        f='r4'
        c='gray'
        r='-r 0,1'
    elif ftype == 'mag':
        q='mag'
        f='cr4'
        c='gray'
    else:
        print('wrong ftype - choose one of: unw,pha,coh,mag')
        return
    outfile = binfile+'.png'
    runcmd('cpxfiddle -w {0} -o sunraster -q {1} -f {2} -c {3} {4} {5} > {6}'.format(str(width), q, f, c, r, binfile, binfile+'.ras'))
    runcmd('convert -resize 700x {0} {1}'.format(binfile+'.ras',outfile))
    return outfile


def main_prepare_masks(ifg, coh, landmask, cohthres, outmask_inside, outmask_full):
    #print('preparing masks (coh-based+landmask, outarea)')
    print('preparing masks (coh-based, landmask+outarea)')
    prefix = os.path.dirname(outmask_inside)+'/'
    tmpmask1 = prefix+'tmp.landtmp.nc'
    tmpmask2 = prefix+'tmp.outtmp.nc'
    #tmpmask3 = prefix+'tmp.intmp.nc'
    #cut landmask to the ifg extents
    runcmd('gmt grdcut -N1 {0} -G{1} -R{2}; gmt grdedit {1} -T -R{2}'.format(landmask, tmpmask1, ifg))
    runcmd('gmt grdmath {0} 0 NAN 0 GT 0 DENAN = {1}'.format(coh, tmpmask2))
    # runcmd('gmt grdmath {0} 0 NAN {1} GT 1 DENAN = {2}'.format(coh, str(cohthres), tmpmask3))
    #ok, use only coh based masking..
    runcmd('gmt grdmath {0} 0 NAN {1} GT 1 DENAN = {2}'.format(coh, str(cohthres), outmask_inside))
    #runcmd('gmt grdmath -N {0} {1} MUL = {2}'.format(tmpmask3, tmpmask1, outmask_inside))
    runcmd('gmt grdmath -N {0} {1} MUL {2} MUL = {3}'.format(tmpmask2, tmpmask1, outmask_inside, outmask_full))
    #runcmd('gmt grdmath -N {0} {1} MUL = {2}'.format(tmpmask2, outmask_inside, outmask_full))
    #cleaning
    os.remove(tmpmask1)
    os.remove(tmpmask2)


def resize_bin(inbin, inwid, inlen, outbin, outwid, outlen, dtype = np.byte, intertype = cv2.INTER_NEAREST):
    a = np.fromfile(inbin, dtype=dtype).reshape(inlen, inwid)
    #use cv2.INTER_CUBIC for upsample np.float32 data ...
    #a = a - np.nanmedian(a)
    out = cv2.resize(a,dsize=(outwid,outlen), interpolation=intertype)
    out.astype(dtype).tofile(outbin)
    #os.remove(tmpmask3)


def main_gapfill_ifg(ifg, mask_in, outifgfilled):
    print('gapfilling masked areas')
    prefix = os.path.dirname(outifgfilled)+'/'
    #convert mask to geo grid registration
    tmpmask = prefix+'tmp.mask.nc'
    tmptofill = prefix+'tmp.tofill.nc'
    tmpout = prefix+'tmp.out.nc'
    runcmd('cp {0} {1}; gmt grdedit {1} -T -R{2}'.format(mask_in, tmpmask, ifg))
    runcmd('gmt grdmath -N {0} 0 NAN 10 DENAN {1} MUL 0 NAN = {2}'.format(ifg, tmpmask, tmptofill))
    #filling itself
    runcmd('gmt grdfill {0} -An -G{1}'.format(tmptofill, tmpout))
    runcmd('gmt grdmath {0} 10 NAN 0 NAN 0 DENAN = {1}'.format(tmpout, outifgfilled))
    os.remove(tmpmask)
    os.remove(tmptofill)
    os.remove(tmpout)


def main_prepare_for_unwrapping(pha, coh, maskin, maskfullin, cpx):
    # R=mag*cos(pha), I=mag*sin(pha)  .... use coh as mag
    #first divide coh files by 255
    prefix = os.path.dirname(cpx)+'/'
    tmpcoh = prefix+'tmp.coh.nc'
    tmppha = prefix+'tmp.pha.nc'
    tmpreal = prefix+'R.nc'
    tmpimag = prefix+'I.nc'
    binR = prefix+'R.bin'
    binI = prefix+'I.bin'
    binCPX = prefix+'cpxifg.bin'
    bincoh = prefix+'coh.bin'
    binmask = prefix+'mask.bin'
    binmaskfull = prefix+'mask.full.bin'
    #
    runcmd('gmt grdmath {0} 255 DIV = {1}'.format(coh, tmpcoh)) #out in px reg
    runcmd('cp {0} {1}; gmt grdedit {1} -T -R{2}'.format(pha, tmppha, tmpcoh))  #change the pha to px reg if it is in gridline (it was!)
    # now finally convert to... R,I:
    runcmd('gmt grdmath {0} COS {1} MUL = {2}'.format(tmppha, tmpcoh, tmpreal))
    runcmd('gmt grdmath {0} SIN {1} MUL = {2}'.format(tmppha, tmpcoh, tmpimag))
    runcmd('gmt grd2xyz -ZTLf -bof {0} > {1}'.format(tmpreal, binR))
    runcmd('gmt grd2xyz -ZTLf -bof {0} > {1}'.format(tmpimag, binI))
    print('preparing data for snaphu')
    # and to cpx binary
    RI2cpx(binR, binI, binCPX)
    # ok, let's also add coh and mask for snaphu
    runcmd('gmt grd2xyz -ZTLf -bof {0} > {1}'.format(tmpcoh, bincoh))
    #and export mask to 0,1 binary - using only the landmask now..
    runcmd('gmt grd2xyz -ZTLc -bof {0} > {1}'.format(maskin, binmask))
    runcmd('gmt grd2xyz -ZTLc -bof {0} > {1}'.format(maskfullin, binmaskfull))
    # (and this is good as we can filter R, I using Gaussian or other filter prior to continuation)
    #os.remove(tmpcoh)


def RI2cpx(R, I, cpxfile):
    r = np.fromfile(R, dtype=np.float32)
    i = np.fromfile(I, dtype=np.float32)
    cpx = np.zeros(len(r)+len(i))
    cpx[0::2] = r
    cpx[1::2] = i
    cpx.astype(np.float32).tofile(cpxfile)


def multilook_bin(binfile, outbinfile, width, mlfactor, dtype = 'cr4'):
    runcmd('cpxfiddle -w {0} -o float -q normal -f {1} -M {2}/{2} {3} > {4}'.format(str(width), dtype, str(mlfactor), binfile, outbinfile))


def coh_from_cpxbin(cpxbin, cohbin, width):
    runcmd('cpxfiddle -w {0} -o float -q mag -f cr4 {1} > {2}'.format(str(width), cpxbin, cohbin))



def main_unwrap(cpxbin, cohbin, maskbin = None, outunwbin = 'unwrapped.bin', width = 0, est = None, bin_pre_remove = None, defomax = 0.6, printout = True):
    if width == 0:
        print('error - width is zero')
        return False
    if bin_pre_remove and est:
        # we will remove phase from the est, prior to processing
        est_np = np.fromfile(est, dtype=np.float32)
        est_rem = np.fromfile(bin_pre_remove, dtype=np.float32)
        est_np = est_np - est_rem
        est_np.tofile(est)
        est_np = ''
        est_rem = ''
    if printout:
        print('processing by snaphu')
    prefix = os.path.dirname(cpxbin)+'/'
    snaphuconffile = make_snaphu_conf(prefix, defomax)
    extracmd = ''
    if est:
        extracmd = "-e {}".format(est)
    starttime = time.time()
    if not maskbin:
        snaphucmd = 'snaphu -f {0} -o {1} -c {2} {3} {4} {5}'.format(snaphuconffile, outunwbin, cohbin, extracmd, cpxbin, str(width))
    else:
        snaphucmd = 'snaphu -f {0} -M {1} -o {2} -c {3} {4} {5} {6}'.format(snaphuconffile, maskbin, outunwbin, cohbin, extracmd, cpxbin, str(width))
    runcmd(snaphucmd, printout) #True)
    if printout:
        elapsed_time = time.time()-starttime
        hour = int(elapsed_time/3600)
        minite = int(np.mod((elapsed_time/60),60))
        sec = int(np.mod(elapsed_time,60))
        print("\nElapsed time: {0:02}h {1:02}m {2:02}s".format(hour,minite,sec))
    return


def create_preview(infile, ftype = 'unwrapped'):
    if 'wrapped' not in ftype:
        print('wrong ftype')
        return False
    if ftype == 'wrapped':
        extra = '0'
    else:
        extra = ''
    tosource = os.path.join(os.environ['LiCSARpath'],'lib','LiCSAR_bash_lib.sh')
    command = 'create_preview_'+ftype
    os.system('source {0}; {1} {2} {3}'.format(tosource, command, infile, extra))


def bin2nc(binfile, masknc, outnc, dtype = np.float32):
    #this will use mask nc to both mask and adapt lat,lon
    unw1 = np.fromfile(binfile,dtype=dtype)
    mask1 = xr.open_dataset(masknc)
    a = mask1.copy(deep=True)
    unw1 = unw1.reshape(a.z.shape)
    unw1 = np.flip(unw1,axis=0)
    #
    a.z.values = unw1*mask1.z.values
    a.z.values[a.z.values==0] = np.nan
    a.z.values = a.z.values - np.nanmedian(a.z.values)
    a.to_netcdf(outnc)


def make_snaphu_conf(sdir, defomax = 1.2):
    snaphuconf = ('STATCOSTMODE  DEFO\n',
            'INFILEFORMAT  COMPLEX_DATA\n',
            'CORRFILEFORMAT  FLOAT_DATA\n'
            'OUTFILEFORMAT FLOAT_DATA\n',
            'ESTFILEFORMAT FLOAT_DATA\n',
            'DEFOMAX_CYCLE '+str(defomax)+'\n',
            'RMTMPTILE TRUE\n')
    snaphuconffile = os.path.join(sdir,'snaphu.conf')
    with open(snaphuconffile,'w') as f:
        for l in snaphuconf:
            f.write(l)
    return snaphuconffile


def make_gacos_ifg(frame, pair, outnc):
    pubdir = os.environ['LiCSAR_public']
    geoframedir = os.path.join(pubdir,str(int(frame[:3])),frame)
    for epoch in pair.split('_'):
        epochgacos = os.path.join(geoframedir,'epochs',epoch,epoch+'.sltd.geo.tif')
        if not os.path.exists(epochgacos):
            return False
    epoch1 = pair.split('_')[0]
    epoch2 = pair.split('_')[1]
    gacos1 = os.path.join(geoframedir,'epochs',epoch1,epoch1+'.sltd.geo.tif')
    gacos2 = os.path.join(geoframedir,'epochs',epoch2,epoch2+'.sltd.geo.tif')
    cmd = 'gmt grdmath {0} {1} SUB = {2}'.format(gacos2, gacos1, outnc)
    print(cmd)
    rc = os.system(cmd)
    if os.path.exists(outnc):
        return outnc
    else:
        print('error in GACOS processing of pair '+pair)
        return False


def remove_height_corr(ifg_ml, corr_thres = 0.5, tmpdir = os.getcwd(), dounw = False):
    ifg_mlc = ifg_ml.copy(deep=True)
    ifg_mlc['toremove'] = 0*ifg_mlc['coh']
    thisisit, thistype = correct_hgt(ifg_mlc, blocklen = 20, tmpdir = tmpdir, dounw = dounw)
    #thisisit can be either False, xr.DataArray, or np.float - ok, adding 'thistype' that can be bool, float, xr
    if not thistype == 'bool':
        if thistype == 'float':
            print('we use average value of {} rad/m'.format(str(thisisit)))
        else:
            print('using hgt correlation grid to reduce hgt component')
        ifg_mlc['toremove'].values = thisisit*ifg_mlc['hgt']
    # first, correlate ifg_ml['pha'] and ifg_ml['hgt'] in blocks - perhaps good to partially unwrap it here, per block - so... snaphu in python!?!
    # get coefficient for correction of correlating areas
    # interpolate the coefficient throughout whole raster
    # multiply by hgt = 'to_remove'
    return ifg_mlc['toremove']


def block_hgtcorr(cpx, coh, hgt, procdir = os.getcwd(), dounw = True, block_id=None):
    #print(block_id)
    #return np.array([[0]])
    # first unwrap the block if conditions are ok
    toret = None
    if hgt[hgt!=0].size == 0:
        toret = np.nan
    if not toret:
        if np.max(hgt[hgt!=0]) - np.min(hgt[hgt!=0]) < 100:
            toret = np.nan
    if not toret:
        # too small coherence... although it may work anyway..
        if np.max(coh)<0.1:
            toret = np.nan
    #if not toret:
    #    if np.mean(coh)<0.05:
    #        toret = np.nan
    if not toret:
        if not dounw:
        #ok, let's do it without the unwrapping first....
            try:
                x = np.ravel(hgt)
                y = np.ravel(np.angle(cpx))
                y = y[x!=0]
                x = x[x!=0]
                if x.size == 0:
                    toret = np.nan
                else:
                    if abs(np.corrcoef(x, y)[0,1]) > 0.3:
                        huber = HuberRegressor()
                        rc = huber.fit(x.reshape(-1,1), y)
                        slope = huber.coef_[0]
                        #slope = np.polyfit(np.ravel(hgt), np.ravel(np.angle(cpx)), deg=1)[0]
                        toret = slope
                    else:
                        toret = np.nan
            except:
                print('error during pha x hgt corr')
                toret = np.nan
        else:
            tmpdir = os.path.join(procdir, str(block_id[0])+'_'+str(block_id[1]))
            if not os.path.exists(tmpdir):
                os.mkdir(tmpdir)
            unwr = unwrap_np(cpx, coh, defomax = 0.6, tmpdir=tmpdir)
            # ok, then correlate - and if higher then 0.4, do linear regression to get rad/m
            try:
                x = np.ravel(hgt)
                y = np.ravel(unwr)
                y = y[x!=0]
                x = x[x!=0]
                if x.size == 0:
                    toret = np.nan
                else:
                    if abs(np.corrcoef(x, y)[0,1]) > 0.4:
                        huber = HuberRegressor()
                        rc = huber.fit(x.reshape(-1,1), y)
                        slope = huber.coef_[0]
                        #slope = np.polyfit(x, y, deg=1)[0]
                        toret = slope
                    else:
                        toret = np.nan
            except:
                print('error during unw x hgt corr')
                toret = np.nan
    return np.array([[toret]])
    #print(tmpdir)
    #out = unw
    #return out


def unwrap_np(cpx, coh, defomax = 0.6, tmpdir=os.getcwd()):
    '''
    unwraps given arrays
    '''
    bincoh = os.path.join(tmpdir,'coh.bin')
    binR = os.path.join(tmpdir,'R.bin')
    binI = os.path.join(tmpdir,'I.bin')
    binCPX = os.path.join(tmpdir,'cpxifg.bin')
    unwbin = os.path.join(tmpdir,'unw.bin')
    # create R, I -> CPX as expected by snaphu
    r = np.real(cpx).tofile(binR)
    i = np.imag(cpx).tofile(binI)
    RI2cpx(binR, binI, binCPX)
    #and coh
    coh.tofile(bincoh)
    # unwrap it
    width = coh.shape[0]
    #with nostdout():
    main_unwrap(binCPX, bincoh, outunwbin = unwbin, width = width, defomax = defomax, printout = False)
    # and load it back
    unw1 = np.fromfile(unwbin,dtype=np.float32)
    unw1 = unw1.reshape(coh.shape)
    #shutil - delete tmpdir!!!!!
    shutil.rmtree(tmpdir)
    return unw1


import dask.array as da

def correct_hgt(ifg_ml, blocklen = 20, tmpdir = os.getcwd(), dounw = True, num_workers = 1):
    #ifg_ml['hgtcorr'] = ifg_ml['pha']
    winsize = (blocklen, blocklen)
    cohb = da.from_array(ifg_ml['coh'].where(ifg_ml.gauss_coh>0.4).astype(np.float32).fillna(0.001), chunks=winsize)
    #phab = da.from_array(ifg_ml['pha'].astype(np.float32).fillna(0), chunks=winsize)
    cpxb = da.from_array(ifg_ml['cpx'].where(ifg_ml.gauss_coh>0.4).astype(np.complex64).fillna(0), chunks=winsize)
    hgtb = da.from_array(ifg_ml['hgt'].where(ifg_ml.gauss_coh>0.4).astype(np.float32).fillna(0), chunks=winsize)
    f = da.map_blocks(block_hgtcorr, cpxb, cohb, hgtb, procdir = tmpdir, dounw = dounw, meta=np.array(()), chunks = (1,1))
    try:
        #with nostdout():
        hgtcorr =  f.compute(num_workers=num_workers)
    except:
        print('error in computing hgt correlation grid')
        return False
    # make it to xr:
    aaa = xr.Dataset()
    aaa['coh'] = ifg_ml['coh'].fillna(0).coarsen({'lat': blocklen, 'lon': blocklen}, boundary='pad').mean()
    aaa['hgtcorr'] = aaa['coh']
    aaa.hgtcorr.values = hgtcorr
    #count means - number of non-nan data..!
    if aaa.hgtcorr.count() > 50:
        kernel = Gaussian2DKernel(x_stddev=1)
        #interpolate nans using gaussian 2d kernel.. lower stddev, and iterate till all nans are replaced!
        while aaa.hgtcorr.isnull().max() == True:
            aaa.hgtcorr.values = interpolate_replace_nans(aaa.hgtcorr.values, kernel)
        #interpolate it to the higher resolution
        out = aaa.hgtcorr.interp_like(ifg_ml, method='linear')
        # but here edges are again nans!
        kernel = Gaussian2DKernel(x_stddev=1.5)
        while out.isnull().max() == True:
            out.values = interpolate_replace_nans(out.values, kernel)
        outype = 'xr'
    else:
        out = np.nanmedian(hgtcorr)
        outype = 'float'
        if np.isnan(out):
            print('all NaNs in hgt corr')
            return False, 'bool'
        if np.abs(out) < 0.001:
            print('almost nothing to reduce for hgt')
            return False, 'bool'
    #if we got here, means, now it is up to splining it to the full (ml) resolution... but --- at this moment, i will just use average value
    return out, outype


def process_ifg(frame, pair, 
        procdir = os.getcwd(), 
        ml = 10, fillby = 'gauss', 
        thres = 0.35, prevest = None, 
        remove_pha = None, hgtcorr = False, pre_detrend=True,
        gacoscorr = True, outtif = None, 
        defomax = 0.6, add_resid = True):
    pubdir = os.environ['LiCSAR_public']
    geoframedir = os.path.join(pubdir,str(int(frame[:3])),frame)
    geoifgdir = os.path.join(geoframedir,'interferograms',pair)
    tmpdir = os.path.join(procdir,pair,'temp_'+str(ml))
    tmpgendir = os.path.join(procdir,pair,'temp_gen')
    if not os.path.exists(os.path.join(procdir,pair)):
        os.mkdir(os.path.join(procdir,pair))
    if not os.path.exists(tmpdir):
        os.mkdir(tmpdir)
    if not os.path.exists(tmpgendir):
        os.mkdir(tmpgendir)
    #orig files
    # will use only the filtered ifgs now..
    ifg_pha_file = os.path.join(geoifgdir,pair+'.geo.diff_pha.tif')
    coh_file = os.path.join(geoifgdir,pair+'.geo.cc.tif')
    landmask_file = os.path.join(geoframedir,'metadata',frame+'.geo.landmask.tif')
    hgtfile = os.path.join(geoframedir,'metadata',frame+'.geo.hgt.tif')
    #orig file ncs - I do this only because i was too lazy in importing from rio or gdal directly! need to change this!
    origphanc = os.path.join(tmpgendir,'origpha.nc')
    origcohnc = os.path.join(tmpgendir,'origcoh.nc')
    landmasknc = os.path.join(tmpgendir,'landmask.nc')
    hgtnc = os.path.join(tmpgendir,'hgt.nc')
    # do gacos if exists
    if gacoscorr:
        gacoscorrfile = os.path.join(tmpgendir,'gacos.nc')
        try:
            gacoscorrfile = make_gacos_ifg(frame, pair, gacoscorrfile)
        except:
            print('error processing gacos data for pair '+pair)
            gacoscorrfile = False
    else:
        gacoscorrfile = False
    if gacoscorrfile:
        print('GACOS data found, using to improve unwrapping')
        ingacos = xr.open_dataset(gacoscorrfile)
    else:
        gacoscorr = False
    # prepare input dataset and load
    if not os.path.exists(origphanc):
        convert_tif2nc(ifg_pha_file,origphanc)
    if hgtcorr:
        if not os.path.exists(hgtnc):
            convert_tif2nc(hgtfile,hgtnc)
    if not os.path.exists(landmasknc):
        convert_tif2nc(landmask_file,landmasknc)
    if not os.path.exists(origcohnc):
        convert_tif2nc(coh_file,origcohnc, coh = True)
    #use orig ifg and coh
    inpha = xr.open_dataset(origphanc)
    incoh = xr.open_dataset(origcohnc)
    inmask = incoh.copy(deep=True)
    #inhgt = incoh.copy(deep=True)
    inmask.z.values = np.byte(incoh.z > 0)
    if os.path.exists(landmasknc):
        landmask = xr.open_dataset(landmasknc)
        inmask.z.values = landmask.z.values * inmask.z.values
    else:
        landmask = None
    # to get all in on encapsulement (with coords..)
    ifg = inpha.copy(deep=True)
    ifg = ifg.rename({'z':'pha'})
    ifg['coh'] = ifg.pha
    ifg['mask'] = ifg.pha
    #ifg['origpha'] = ifg.pha
    ifg['coh'].values = incoh.z.values
    ifg['mask'].values = inmask.z.values
    if hgtcorr:
        if os.path.exists(hgtnc):
            inhgt = xr.open_dataset(hgtnc)
        ifg['hgt'] = ifg.pha
        try:
            ifg['hgt'].values = inhgt.z.values
        except:
            print('ERROR in importing heights!')
    if gacoscorr:
        ifg['gacos'] = ifg.pha
        ifg['gacos'].values = ingacos.z.values
    # erase from memory
    inpha = ''
    incoh = ''
    inmask = ''
    ingacos = ''
    inhgt = ''
    # now doing multilooking, using coh as mag...
    #make complex from coh and pha
    ifg['cpx'] = ifg.coh.copy()
    ifg['cpx'].values = magpha2RI_array(ifg.coh.values, ifg.pha.values)
    #fixing difference in xarray version... perhaps...
    if 'lat' not in ifg.coords:
        print('warning - perhaps old xarray version - trying anyway')
        ifg = ifg.rename_dims({'x':'lon','y':'lat'})
    #ifg_ml = multilook_by_gauss(ifg, ml)
    # result more or less same, but much more economic way in memory..
    #WARNING - ONLY THIS FUNCTION HAS GACOS INCLUDED NOW! (and heights fix!!!)
    ifg_ml = multilook_normalised(ifg, ml, tmpdir = tmpdir, hgtcorr = hgtcorr, pre_detrend = pre_detrend)
    #here we should keep the (not wrapped) phase we remove due to corrections - here, heights, and gacos
    #mask it ... this worked ok here:
    #thres = 0.5
    #gmmm.... it works better if i use the 'trick on gauss of normalised ifg....'
    ifg_ml['mask_coh'] = ifg_ml['mask'].where(ifg_ml.gauss_coh > thres).fillna(0)
    ifg_ml['origpha'] = ifg_ml.pha.copy(deep=True)
    #try without modal filter..
    #ifg_ml = filter_mask_modal(ifg_ml, 'mask_coh', 'mask_coh', 8)
    print('interpolate coh-based masked areas of gauss pha')
    #tofillpha = ifg_ml.pha.fillna(0).where(ifg_ml.mask_coh.where(ifg_ml.mask == 1).fillna(1) == 1)
    if fillby == 'gauss':
        #trying astropy approach now:
        print('filling through Gaussian kernel')
        kernel = Gaussian2DKernel(x_stddev=2)
        # create a "fixed" image with NaNs replaced by interpolated values
        tempar_mag1 = np.ones_like(ifg_ml.pha)
        cpxarr = magpha2RI_array(tempar_mag1, ifg_ml.pha.fillna(0).values)
        ifg_ml['cpx_tofill'] = ifg['cpx']
        ifg_ml['cpx_tofill'].values = cpxarr
        tofill = ifg_ml['cpx_tofill'].where(ifg_ml.mask_coh.where(ifg_ml.mask == 1).fillna(1) == 1)
        tofillR = np.real(tofill)
        tofillI = np.imag(tofill)
        filledR = interpolate_replace_nans(tofillR.values, kernel)
        filledI = interpolate_replace_nans(tofillI.values, kernel)
        ifg_ml['pha'].values = np.angle(filledR + 1j*filledI)
        #sometimes the whole area is not within gauss kernel - use NN for that:
        ifg_ml['pha'].values = interpolate_nans(ifg_ml['pha'].values, method='nearest')
        #ifg_ml.pha.fillna(0).where(ifg_ml.mask_coh.where(ifg_ml.mask == 1).fillna(1) == 1)
        #ifg_ml['pha'].values = interpolate_replace_nans(tofillpha.values, kernel)
    else:
        print('filling by nearest neighbours')
        #tofillpha = ifg_ml.pha.fillna(0).where(ifg_ml.mask_coh.where(ifg_ml.mask == 1).fillna(1) == 1)
        tofillpha = ifg_ml.pha.fillna(0).where(ifg_ml.mask_coh.where(ifg_ml.mask == 1).fillna(1) == 1)
        ifg_ml['pha'].values = interpolate_nans(tofillpha.values, method='nearest')
        #ifg_ml['gauss_pha'] = ifg_ml['gauss_pha'].fillna(0)
    #exporting for snaphu
    #normalise mag from the final pha
    tempar_mag1 = np.ones_like(ifg_ml.pha)
    cpxarr = magpha2RI_array(tempar_mag1, ifg_ml.pha.fillna(0).values) #no need to fillna, but just in case...
    ifg_ml['gauss_cpx'].values = cpxarr
    print('unwrapping by snaphu')
    binmask= os.path.join(tmpdir,'gaussmask.bin')
    #bincoh = os.path.join(tmpdir,'gausscoh.bin')
    bincoh = os.path.join(tmpdir,'coh.bin')
    binR = os.path.join(tmpdir,'gaussR.bin')
    binI = os.path.join(tmpdir,'gaussI.bin')
    binCPX = os.path.join(tmpdir,'cpxgaussifg.bin')
    outunwbin = os.path.join(tmpdir,'gaussunwrapped.bin')
    #print('exporting to bin files')
    ifg_ml.mask_coh.fillna(0).values.astype(np.byte).tofile(binmask)
    #ifg_ml.gauss_coh.fillna(0.001).values.astype(np.float32).tofile(bincoh)
    # we should use the orig coh for weights... and perhaps very low coh values instead of 0
    ifg_ml.coh.fillna(0.001).values.astype(np.float32).tofile(bincoh)
    r = np.real(ifg_ml.gauss_cpx).astype(np.float32).fillna(0).values.tofile(binR)
    i = np.imag(ifg_ml.gauss_cpx).astype(np.float32).fillna(0).values.tofile(binI)
    RI2cpx(binR, binI, binCPX)
    width = len(ifg_ml.lon)
    length = len(ifg_ml.lat)
    #unwrapping itself
    if type(prevest) != type(None):
        #resizing previous ML step and using to unwrap
        bin_pre = os.path.join(tmpdir,'prevest.bin')
        bin_est = os.path.join(tmpdir,'prevest.rescaled.bin')
        bin_pre_remove = os.path.join(tmpdir,'prevest.rescaled.remove.bin')
        #binI_pre = os.path.join(tmpdir,'prevest.I.bin')
        
        kernel = Gaussian2DKernel(x_stddev=2)
        prevest.values = interpolate_replace_nans(prevest.where(prevest != 0).values, kernel)
        #prevest.values = interpolate_replace_nans(prevest.values, kernel)
        #filling other nans to 0 - this can happen if null regions are too large (larger than the kernel accepts)
        prevest.astype(np.float32).fillna(0).values.tofile(bin_pre)
        width_pre = len(prevest.lon)
        length_pre = len(prevest.lat)
        #resize_bin(bin_pre, width_pre, length_pre, bin_est, width, length, dtype = np.float32, intertype = cv2.INTER_CUBIC)
        resize_bin(bin_pre, width_pre, length_pre, bin_est, width, length, dtype = np.float32, intertype = cv2.INTER_LINEAR)
        ifg_ml['toremove'].astype(np.float32).fillna(0).values.tofile(bin_pre_remove)
        main_unwrap(binCPX, bincoh, binmask, outunwbin, width, bin_est, bin_pre_remove = bin_pre_remove, defomax = defomax)
    else:
        main_unwrap(binCPX, bincoh, binmask, outunwbin, width, defomax = defomax)
    print('importing snaphu result to ifg_ml')
    binfile = outunwbin
    #toxr = ifg_ml
    daname = 'unw'
    dtype = np.float32
    unw1 = np.fromfile(binfile,dtype=dtype)
    unw1 = unw1.reshape(ifg_ml.pha.shape)
    #unw1 = np.flip(unw1,axis=0)
    ifg_ml[daname] = ifg_ml['pha'].copy(deep=True)
    ifg_ml[daname].values = unw1
    #ok, so the gauss-based coh mask is not the best to do... so exporting 'all pixels'
    ifg_ml[daname] = ifg_ml[daname]*ifg_ml['mask']
    #print('20210722 - testing now - using gauss-based coh mask, ignore the next message:')
    #ifg_ml[daname] = ifg_ml[daname]*ifg_ml['mask_coh']
    ifg_ml[daname].values = ifg_ml[daname].values - np.nanmedian(ifg_ml[daname].where(ifg_ml.mask_coh>0).values)
    #ifg_ml[daname] = ifg_ml[daname]*ifg_ml.mask_coh
    ifg_ml[daname] = ifg_ml[daname]*ifg_ml.mask
    #print('unwrap also residuals from the filtered cpx, and add to the final unw - mask only waters..')
    cpxarr = magpha2RI_array(tempar_mag1, ifg_ml.origpha.fillna(0).where(ifg_ml.mask == 1).values)
    ifg_ml['origcpx'] = ifg_ml['gauss_cpx'].copy(deep=True)
    ifg_ml['origcpx'].values = cpxarr
    if add_resid:
        print('unwrapping residuals and adding back to the final unw output')
        ifg_ml['resid_cpx'] = ifg_ml.origcpx * ifg_ml.gauss_cpx.conjugate()
        incpx = 'resid_cpx'
        binR = os.path.join(tmpdir,incpx+'.R.bin')
        binI = os.path.join(tmpdir,incpx+'.I.bin')
        binCPX = os.path.join(tmpdir,incpx+'.cpx.bin')
        outunwbin = os.path.join(tmpdir,incpx+'.unw.bin')
        binfile = outunwbin
        daname = incpx+'.unw'
        #
        r = np.real(ifg_ml[incpx]).astype(np.float32).fillna(0).values.tofile(binR)
        i = np.imag(ifg_ml[incpx]).astype(np.float32).fillna(0).values.tofile(binI)
        RI2cpx(binR, binI, binCPX)
        main_unwrap(binCPX, bincoh, binmask, outunwbin, width, defomax = 0.5)
        unw1 = np.fromfile(binfile,dtype=dtype)
        unw1 = unw1.reshape(ifg_ml.pha.shape)
        ifg_ml[daname] = ifg_ml['pha'].copy()
        ifg_ml[daname].values = unw1
        ifg_ml[daname] = ifg_ml[daname]*ifg_ml['mask']
        print('debug - avoiding median correction now')
        #ifg_ml[daname].values = ifg_ml[daname].values - np.nanmedian(ifg_ml[daname].where(ifg_ml.mask_coh>0).values)
        #print('again, masking the final product by gauss coh threshold')
        #ifg_ml[daname] = ifg_ml[daname]*ifg_ml.mask_coh
        #
        #
        #ifg_ml['resid_'].plot()
        #now the unw will have the residual phase added back
        ifg_ml['unw'] = ifg_ml['unw'] + ifg_ml[daname]
    #if gacoscorr:
    #    #we have removed GACOS estimate from unw, now time to add it back!
    #    ifg_ml['unw'] = ifg_ml['unw'] + ifg_ml['gacos']
    # add back what we have removed before..
    ifg_ml['unw'] = ifg_ml['unw'] + ifg_ml['toremove']
    #mask it
    ifg_ml['unw'] = ifg_ml['unw'].where(ifg_ml.mask>0)
    # ok ok.... let's mask by the gauss mask.. although, can be quite missing lot of areas
    #ifg_ml['unw'] = ifg_ml['unw'].where(ifg_ml.mask_coh>0)
    # hmmm... just to make it nicer...
    ifg_ml['unw'] = ifg_ml['unw'] - ifg_ml['unw'].median()
    if outtif:
        #ifg_ml.pha.fillna(0).where(ifg_ml.mask_coh.where(ifg_ml.mask == 1).fillna(1) == 1)
        #ifg_ml['unw'].to_netcdf(outtif+'.nc')
        ifg_ml['unw'].to_netcdf(outtif+'.nc')
        rc = os.system('gmt grdconvert -G{0}=gd:GTiff -R{1} {0}.nc'.format(outtif, ifg_pha_file))
        #rc = os.system('source {0}/lib/LiCSAR_bash_lib.sh; create_preview_unwrapped {1} {2}'.format(os.environ['LiCSARpath'], outtif, frame))
        rc = os.system('source {0}/lib/LiCSAR_bash_lib.sh; create_preview_unwrapped {1}'.format(os.environ['LiCSARpath'], outtif))
        try:
            os.remove(outtif+'.nc')
        except:
            print('ERROR removing the nc file - something wrong with export')
    return ifg_ml


def process_frame(frame, ml = 10, hgtcorr = False, cascade=False):
    if cascade and ml>1:
        print('error - the cascade approach is ready only for ML1')
        return False
    #the best to run in directory named by the frame id
    pubdir = os.environ['LiCSAR_public']
    geoframedir = os.path.join(pubdir,str(int(frame[:3])),frame)
    geoifgdir = os.path.join(geoframedir,'interferograms')
    hgtfile = os.path.join(geoframedir,'metadata', frame+'.geo.hgt.tif')
    raster = gdal.Open(hgtfile)
    framewid = raster.RasterXSize
    framelen = raster.RasterYSize
    for pair in os.listdir(geoifgdir):
        if os.path.exists(os.path.join(geoifgdir, pair, pair+'.geo.diff_pha.tif')):
            #check its dimensions..
            raster = gdal.Open(os.path.join(geoifgdir, pair, pair+'.geo.diff_pha.tif'))
            if (framewid != raster.RasterXSize) or (framelen != raster.RasterYSize):
                #use tolerance of max pixels
                maxpixels = 4
                if ((abs(framewid - raster.RasterXSize) > maxpixels) or (abs(framelen - raster.RasterYSize) > maxpixels)):
                    print('ERROR - the file {} has unexpected dimensions, skipping'.format(os.path.join(geoifgdir, pair, pair+'.geo.diff_pha.tif')))
                    continue
                print('ERROR - the file {} has unexpected dimensions, trying to fix'.format(os.path.join(geoifgdir, pair, pair+'.geo.diff_pha.tif')))
                for tif in glob.glob(os.path.join(geoifgdir, pair, pair+'.geo.*.tif')):
                    outfile = tif+'.tmp.tif'
                    try:
                        filedone = reproject_to_match(tif, hgtfile, outfile)
                        if os.path.exists(outfile):
                            shutil.move(outfile, tif)
                    except:
                        print('something wrong during reprojection, skipping')
                        continue
                    #os.system('gmt grdsample {0} -G{1}')
            try:
                raster = gdal.Open(os.path.join(geoifgdir, pair, pair+'.geo.diff_pha.tif'))
                if (framewid != raster.RasterXSize) or (framelen != raster.RasterYSize):
                    print('ERROR - the file {} has unexpected dimensions, skipping'.format(os.path.join(geoifgdir, pair, pair+'.geo.diff_pha.tif')))
                    continue
            except:
                print('some error processing file {}'.format(os.path.join(geoifgdir, pair, pair+'.geo.diff_pha.tif')))
                continue
            if not os.path.exists(os.path.join(pair,pair+'.unw')):
                print('processing pair '+pair)
                try:
                    if cascade:
                        ifg_ml = cascade_unwrap(frame, pair, procdir = os.getcwd())
                    else:
                        ifg_ml = process_ifg(frame, pair, procdir = os.getcwd(), ml = ml, hgtcorr = hgtcorr, fillby = 'gauss')
                    #else:
                    #    print('ML set to 1 == will try running the cascade (multiscale) unwrapping')
                    #    ifg_ml = cascade_unwrap(frame, pair, procdir = os.getcwd())
                    #ifg_ml.unw.values.tofile(pair+'/'+pair+'.unw')
                    #np.flipud(ifg_ml.unw.fillna(0).values).tofile(pair+'/'+pair+'.unw')
                    #np.flipud(ifg_ml.unw.where(ifg_ml.mask_coh > 0).values).astype(np.float32).tofile(pair+'/'+pair+'.unw')
                    np.flipud(ifg_ml.unw.where(ifg_ml.mask > 0).values).astype(np.float32).tofile(pair+'/'+pair+'.unw')
                    #or use gauss here?
                    np.flipud((ifg_ml.coh.where(ifg_ml.mask > 0)*255).astype(np.byte).fillna(0).values).tofile(pair+'/'+pair+'.cc')
                    width = len(ifg_ml.lon)
                    create_preview_bin(pair+'/'+pair+'.unw', width, ftype = 'unw')
                    os.system('rm '+pair+'/'+pair+'.unw.ras')
                    os.system('rm -r '+pair+'/'+'temp_'+str(ml))
                except:
                    print('ERROR processing of pair '+pair)
                    os.system('rm -r '+pair)
            if not os.path.exists(os.path.join(pair,pair+'.unw')):
                print('some error occured and the unw was not processed')
                os.system('rm -r '+pair)


def multilook_by_gauss(ifg, ml = 10):
    kernel = Gaussian2DKernel(x_stddev=ml)
    resultR = convolve_fft(np.real(ifg['cpx'].values), kernel, allow_huge=True)
    resultI = convolve_fft(np.imag(ifg['cpx'].values), kernel, allow_huge=True)
    ifg['gauss_pha'] = ifg['pha'].copy()
    ifg['gauss_pha'].values = np.angle(resultR + 1j*resultI)
    #only now multilook... using coh as mag
    ifg['gauss_cpx'] = ifg['cpx'].copy()
    ifg['gauss_cpx'].values = magpha2RI_array(ifg['coh'].values, ifg['gauss_pha'].values)
    bagcpx = ifg[['gauss_cpx']].where(ifg.mask>0).coarsen({'lat': ml, 'lon': ml}, boundary='trim')
    ifg_ml = bagcpx.sum() / bagcpx.count()
    #just make it back to pha and coh..
    ifg_ml['pha'] = ifg_ml.gauss_cpx.copy()
    ifg_ml['pha'].values = np.angle(ifg_ml.gauss_cpx)
    #make the coh 'same' as orig (multilooked only)
    ifg_ml['coh'] = ifg_ml['pha'].copy()
    ifg_ml.coh.values = np.abs(ifg_ml.gauss_cpx)
    #downsample mask
    ifg_ml['mask'] = ifg.mask.coarsen({'lat': ml, 'lon': ml}, boundary='trim').max()
    #ifg_ml = make_gauss_coh(ifg_ml)
    ifg_ml['gauss_coh'] = make_gauss_coh_good(ifg, ml = ml)
    return ifg_ml


def make_gauss_coh(ifg_ml):
    tempar_mag1 = np.ones_like(ifg_ml.pha)
    cpxarr = magpha2RI_array(tempar_mag1, ifg_ml.pha.values)
    ifg_ml['cpx'] = ifg_ml.gauss_cpx.copy(deep=True)
    ifg_ml['cpx'].values = cpxarr
    #print('filter using (adapted) gauss filter')
    ifg_ml['temp_cpx'] = filter_cpx_gauss(ifg_ml, sigma = 1.2, trunc = 4)
    #use magnitude after filtering as coherence
    ifg_ml['gauss_coh'] = ifg_ml.pha.copy()
    ifg_ml['gauss_coh'].values = np.abs(ifg_ml.temp_cpx.values)
    ifg_ml = ifg_ml.drop('cpx')
    ifg_ml = ifg_ml.drop('temp_cpx')
    return ifg_ml


def make_gauss_coh_good(ifg, ml = 10):
    ifgtemp = ifg.copy(deep=True)
    ifgtemp['cpx'] = ifgtemp.coh.copy()
    ifgtemp['cpx'].values = magpha2RI_array(ifgtemp.coh.values, ifgtemp.pha.values)
    #landmask it and multilook it
    bagcpx = ifgtemp[['cpx']].where(ifgtemp.mask>0).coarsen({'lat': ml, 'lon': ml}, boundary='trim')
    ifgtemp_ml = bagcpx.sum() / bagcpx.count()
    ifgtemp_ml['pha'] = ifgtemp_ml.cpx.copy()
    ifgtemp_ml['pha'].values = np.angle(ifgtemp_ml.cpx)
    #downsample mask
    ifgtemp_ml['mask'] = ifgtemp.mask.coarsen({'lat': ml, 'lon': ml}, boundary='trim').max()
    #normalise mag
    tempar_mag1 = np.ones_like(ifgtemp_ml.pha)
    cpxarr = magpha2RI_array(tempar_mag1, ifgtemp_ml.pha.values)
    ifgtemp_ml['cpx'].values = cpxarr
    ifgtemp_ml['gauss_cpx'] = filter_cpx_gauss(ifgtemp_ml, sigma = 1, trunc = 2)
    ifgtemp_ml['gauss_pha'] = ifgtemp_ml.pha
    ifgtemp_ml['gauss_pha'].values = np.angle(ifgtemp_ml.gauss_cpx.values)
    #use magnitude after filtering as coherence
    ifgtemp_ml['gauss_coh'] = ifgtemp_ml.pha.copy()
    ifgtemp_ml['gauss_coh'].values = np.abs(ifgtemp_ml.gauss_cpx.values)
    return ifgtemp_ml['gauss_coh']


def multilook_normalised(ifg, ml = 10, tmpdir = os.getcwd(), hgtcorr = True, pre_detrend = True):
    #landmask it and multilook it
    if ml > 1:
        bagcpx = ifg[['cpx']].where(ifg.mask>0).coarsen({'lat': ml, 'lon': ml}, boundary='trim')
        ifg_ml = bagcpx.sum() / bagcpx.count()
    else:
        ifg_ml = ifg[['cpx']].where(ifg.mask>0)
    #prepare 'toremove' layer
    ifg_ml['toremove'] = ifg_ml.cpx.copy()
    ifg_ml['toremove'].values = 0*np.angle(ifg_ml.cpx)
    if pre_detrend:
        ifg_ml['cpx'], ifg_ml['toremove'] = detrend_ifg_xr(ifg_ml['cpx'], isphase=False, return_correction = True)
    ifg_ml['pha'] = ifg_ml.cpx.copy()
    ifg_ml['pha'].values = np.angle(ifg_ml.cpx)
    #have the orig coh here:
    ifg_ml['coh'] = 0*ifg_ml['pha']
    ifg_ml.coh.values = np.abs(ifg_ml.cpx)
    #keep the original original pha values
    ifg_ml['origpha_noremovals'] = ifg_ml.pha.copy(deep=True)
    #downsample mask
    if ml > 1:
        ifg_ml['mask'] = ifg.mask.coarsen({'lat': ml, 'lon': ml}, boundary='trim').max()
        if 'gacos' in ifg.variables:
            ifg_ml['gacos'] = ifg.gacos.coarsen({'lat': ml, 'lon': ml}, boundary='trim').mean()
        if 'hgt' in ifg.variables:
            ifg_ml['hgt'] = ifg.hgt.coarsen({'lat': ml, 'lon': ml}, boundary='trim').mean()
    else:
        ifg_ml['mask'] = ifg.mask
        if 'gacos' in ifg.variables:
            ifg_ml['gacos'] = ifg.gacos
        if 'hgt' in ifg.variables:
            ifg_ml['hgt'] = ifg.hgt
    #have gacos removed first, prior to doing height corr:
    if 'gacos' in ifg_ml.variables:
        ifg_ml['pha'].values = wrap2phase(ifg_ml['pha'] - ifg_ml['gacos'])
        ifg_ml['pha'] = ifg_ml['pha'].where(ifg_ml.mask>0)
        ifg_ml['gacos'] = ifg_ml['gacos'].where(ifg_ml.mask>0)
    ifg_ml['pha'] = ifg_ml['pha'].where(ifg_ml.mask>0)
    #now fix the correlation with heights:
    if hgtcorr:
        ifg_ml = filter_ifg_ml(ifg_ml)
        #ifg_ml['toremove'] = ifg_ml['toremove'] + 
        toremove_hgt = remove_height_corr(ifg_ml, tmpdir = tmpdir)
        ifg_ml['toremove'] = ifg_ml['toremove'] + toremove_hgt
        ifg_ml['pha'].values = wrap2phase(ifg_ml['pha'] - toremove_hgt)
        ifg_ml['pha'] = ifg_ml['pha'].where(ifg_ml.mask>0)
        # ifg_ml['pha'].plot()
    #else:
    #    ifg_ml['toremove'] = 0*ifg_ml['pha']
    #here we remove height correlation - gacos was already removed from pha
    #ifg_ml['pha'].values = wrap2phase(ifg_ml['pha'] - ifg_ml['toremove'])
    #ifg_ml['pha'] = ifg_ml['pha'].where(ifg_ml.mask>0)
    # maybe not the best, but have gacos correction inside the toremove variable
    if 'gacos' in ifg_ml.variables:
        ifg_ml['toremove'] = ifg_ml['toremove'] + ifg_ml['gacos'] #.fillna(0)
    ifg_ml['toremove'] = ifg_ml['toremove'].where(ifg_ml.mask>0)
    #print('filter using (adapted) gauss filter')
    ifg_ml = filter_ifg_ml(ifg_ml)
    return ifg_ml


def load_tif2xr(tif):
    xrpha = xr.open_rasterio(tif)
    xrpha = xrpha.squeeze('band')
    xrpha = xrpha.drop('band')
    xrpha = xrpha.rename({'x': 'lon','y': 'lat'})
    return xrpha


def export_xr2tif(xrda, tif):
    xrda = xrda.astype(np.float32)
    coordsys = xrda.crs.split('=')[1]
    xrda = xrda.rio.set_spatial_dims(x_dim="lon", y_dim="lat", inplace=True)
    xrda = xrda.rio.write_crs(coordsys, inplace=True)
    xrda.rio.to_raster(tif, compress='deflate')
    

def deramp_ifg_tif(phatif, unwrap_after = True):
    # works fine if the path is to some diff_pha.tif file inside $LiCSAR_public only!
    phatif = os.path.realpath(phatif)
    if not os.path.exists(phatif):
        print('the input tif does not exist, exiting')
        return 0
    xrpha = load_tif2xr(phatif)
    xrpha_detrended = detrend_ifg_xr(xrpha, isphase=True)
    xrpha_detrended = xrpha_detrended.where(xrpha != 0).fillna(0)
    phatif_orig = phatif.replace('.geo.diff_pha.tif','.geo.diff_pha.orig.tif')
    rc = os.system('mv {0} {1}'.format(phatif, phatif_orig))
    rc = os.system('mv {0} {1}'.format(phatif.replace('.geo.diff_pha.tif','.geo.diff.png'), 
       phatif.replace('.geo.diff_pha.tif','.geo.diff.orig.png')))
    export_xr2tif(xrpha_detrended, phatif)
    # just do also preview
    create_preview(phatif, 'wrapped')
    if unwrap_after:
        # and probably reunwrapping is needed...
        # doing it original way
        frame = phatif.split('/')[-4]
        pair = phatif.split('/')[-2]
        unwtif = phatif.replace('.geo.diff_pha.tif','.geo.unw.tif')
        unwtif_orig = unwtif.replace('.geo.unw.tif','.geo.unw.orig.tif')
        os.system('mv {0} {1}'.format(unwtif, unwtif_orig))
        os.system('mv {0} {1}'.format(unwtif.replace('.geo.unw.tif','.geo.unw.png'), unwtif.replace('.geo.unw.tif','.geo.unw.orig.png')))
        os.system('unwrap_geo.sh {0} {1}'.format(frame, pair))
        #or using this approach? well... it takes some 5x+ longer!!!, so perhaps not
        #outtif = phatif.replace('.geo.diff_pha.tif','.geo.unw.tif')
        #process_ifg(frame, pair, ml = 1, pre_detrend = False, outtif = outtif)
    
    


def detrend_ifg_xr(xrda, isphase=True, return_correction = False):
    #not done yet... doris matlab cpxdetrend function is too cool.. need to understand it first
    #xrcpx = ifg_ml['cpx']
    da = xrda.copy(deep=True).fillna(0)
    if isphase:
        #convert to cpx values first
        mag = np.ones(da.values.shape)
        cpx = magpha2RI_array(mag, da.values)
    else:
        cpx = da.values
    fftt = np.fft.fft2(cpx)
    fftt = np.abs(fftt)
    #remove zero line and column - often has too much of zeroes there..
    fftt[0] = fftt[0]*0
    fftt = fftt.transpose()
    fftt[0] = fftt[0]*0
    fftt = fftt.transpose()
    numfringesx = np.argmax(np.sum(fftt,axis=0))
    numfringesy = np.argmax(np.sum(fftt,axis=1))
    [Y,X]=fftt.shape
    if numfringesx > X/2:
        numfringesx = numfringesx - X
    if numfringesy > Y/2:
        numfringesy = numfringesy - Y
    print('detrending by {0}/{1} fringes in lon/lat'.format(numfringesx, numfringesy))
    trendx = np.linspace(0,2*np.pi,X) * numfringesx
    trendy = np.linspace(0,2*np.pi,Y) * numfringesy
    trendx = np.tile(trendx, (Y,1))
    trendy = np.tile(trendy, (X,1)).transpose()
    correction = trendx + trendy
    if isphase:
        #was phase, return phase
        da.values = wrap2phase(da.values - correction)
    else:
        mag = np.ones(da.values.shape)
        cpx_corr = magpha2RI_array(mag, correction)
        #a trick to apply the correction only to non-nan values.. 'other'
        #da = da.where(xrda.isnull(), other=da.values * np.conjugate(cpx_corr))
        da.values = da.values * np.conjugate(cpx_corr)
        #da.where(xrcpx.isnull()) = 0
    if return_correction:
        corrda = da.copy()
        corrda.values = correction
        return da, corrda
    else:
        return da



    
def filter_ifg_ml(ifg_ml):
    #normalise mag
    tempar_mag1 = np.ones_like(ifg_ml.pha)
    cpxarr = magpha2RI_array(tempar_mag1, ifg_ml.pha.values)
    ifg_ml['cpx'].values = cpxarr
    print('filter using (adapted) gauss filter')
    ifg_ml['gauss_cpx'] = filter_cpx_gauss(ifg_ml, sigma = 1, trunc = 2)
    ifg_ml['gauss_pha'] = 0*ifg_ml['pha']
    ifg_ml['gauss_pha'].values = np.angle(ifg_ml.gauss_cpx.values)
    #use magnitude after filtering as coherence
    ifg_ml['gauss_coh'] = 0*ifg_ml['pha']
    ifg_ml['gauss_coh'].values = np.abs(ifg_ml.gauss_cpx.values)
    return ifg_ml


def block_unwrap(da_pha, da_coh):
    return da_unw


def block_correlate(da1, da2):
    return dacorr


def wrap2phase(A):
    return np.angle(np.exp(1j*A))


import time

def cascade_unwrap(frame, pair, procdir = os.getcwd(), outtif = None):
    print('performing cascade unwrapping')
    starttime = time.time()
    ifg_ml20 = process_ifg(frame, pair, procdir = procdir, ml = 20, fillby = 'gauss', defomax = 0.5, add_resid = False)
    ifg_ml10 = process_ifg(frame, pair, procdir = procdir, ml = 10, fillby = 'gauss', prevest = ifg_ml20['unw'], defomax = 0.5, add_resid = False)
    #elapsed_time10 = time.time()-starttime
    ifg_ml5 = process_ifg(frame, pair, procdir = procdir, ml = 5, fillby = 'gauss', prevest = ifg_ml10['unw'], defomax = 0.6, add_resid = False)
    ifg_ml3 = process_ifg(frame, pair, procdir = procdir, ml = 3, fillby = 'gauss', prevest = ifg_ml5['unw'], defomax = 0.6, add_resid = False)
    ifg_ml1 = process_ifg(frame, pair, procdir = procdir, ml = 1, fillby = 'gauss', prevest = ifg_ml3['unw'], defomax = 1.0, add_resid = True, outtif=outtif)
    elapsed_time = time.time()-starttime
    hour = int(elapsed_time/3600)
    minite = int(np.mod((elapsed_time/60),60))
    sec = int(np.mod(elapsed_time,60))
    print("\nTotal elapsed time: {0:02}h {1:02}m {2:02}s".format(hour,minite,sec))
    return ifg_ml1

#print('preparing nc and png')
#bin2nc(outunwbin, mask_full, outunwnc, dtype = np.float32)
#create_preview(outunwnc, ftype = 'unwrapped')
