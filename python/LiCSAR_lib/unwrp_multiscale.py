################################################################################
#Imports
################################################################################
import numpy as np
import os, glob
from osgeo import gdal
import subprocess
import xarray as xr
import rioxarray
import pandas as pd

try:
    import cv2
except:
    print('cv2 not loaded - cascade will not work')


from scipy import interpolate
from scipy import ndimage
import time
import matplotlib.pyplot as plt
xr.set_options(keep_attrs=True)


from scipy.ndimage import gaussian_filter
from scipy.ndimage import generic_filter
from scipy import stats
import scipy.signal as sps
import scipy.linalg as spl

from astropy.convolution import Gaussian2DKernel, interpolate_replace_nans, convolve_fft
from sklearn.linear_model import HuberRegressor

import shutil
# fc used only for amp stability..
try:
    import framecare as fc
except:
    print('framecare not loaded')


try:
    from LiCSAR_lib.LiCSAR_misc import *
except:
    print('licsar misc not loaded')


try:
    import LiCSBAS_io_lib as io
    from LiCSBAS_tools_lib import *
except:
    print('licsbas not loaded')

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


def coh_from_phadiff(phadiff, winsize = 3):
    variance = ndimage.generic_filter(phadiff, np.var, size=winsize)
    outcoh = 1/np.sqrt(1+winsize*winsize*variance)
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
    runcmd('cpxfiddle -w {0} -o sunraster -q {1} -f {2} -c {3} {4} {5} > {6} 2>/dev/null'.format(str(width), q, f, c, r, binfile, binfile+'.ras'))
    runcmd('convert -resize 700x {0} {1}'.format(binfile+'.ras',outfile))
    rc = os.system('rm {0}'.format(binfile+'.ras'))
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
    # we may either load R, I from file:
    if type(R) == type('str'):
        r = np.fromfile(R, dtype=np.float32)
        i = np.fromfile(I, dtype=np.float32)
    else:
        r = R.astype(np.float32).ravel()
        i = I.astype(np.float32).ravel()
    cpx = np.zeros(len(r)+len(i))
    cpx[0::2] = r
    cpx[1::2] = i
    cpx.astype(np.float32).tofile(cpxfile)


def multilook_bin(binfile, outbinfile, width, mlfactor, dtype = 'cr4'):
    runcmd('cpxfiddle -w {0} -o float -q normal -f {1} -M {2}/{2} {3} > {4} 2>/dev/null'.format(str(width), dtype, str(mlfactor), binfile, outbinfile))


def coh_from_cpxbin(cpxbin, cohbin, width):
    runcmd('cpxfiddle -w {0} -o float -q mag -f cr4 {1} > {2} 2>/dev/null'.format(str(width), cpxbin, cohbin))


def remove_islands(npa, pixelsno = 50):
    '''
    removes isolated clusters of pixels from numpy array npa having less than pixelsno pixels
    '''
    #check the mask - should be 1 for islands and 0 for nans
    mask = ~np.isnan(npa)
    from scipy import ndimage
    islands, ncomp = ndimage.label(mask)
    for i in range(ncomp):
        #island = islands == i # need to get this one right
        #island = npa[islands==i]
        numofpixels = len(mask[islands==i])
        if numofpixels < pixelsno:
            npa[islands==i] = np.nan
    return npa


def rad2mm_s1(inrad):
    speed_of_light = 299792458 #m/s
    radar_freq = 5.405e9  #for S1
    wavelength = speed_of_light/radar_freq #meter
    coef_r2m = -wavelength/4/np.pi*1000 #rad -> mm, positive is -LOS
    outmm = inrad*coef_r2m
    return outmm


def rad2mm(pha, lam = 0.0312):
    coef_r2m = -lam/4/np.pi*1000 #rad -> mm, positive is -LOS
    outmm = pha*coef_r2m
    return outmm


def plotit():
    plt.figure(figsize=(7, 3.5))
    plt.subplot(1, 2, 1)
    plt.imshow(sig)
    plt.axis('off')
    plt.title('sig')
    plt.subplot(1, 2, 2)
    plt.imshow(mask, cmap=plt.cm.gray)
    plt.axis('off')
    plt.title('mask')
    plt.subplots_adjust(wspace=.05, left=.01, bottom=.01, right=.99, top=.9)
    plt.show()

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


def make_gacos_ifg(frame, pair, outfile):
    print('preparing GACOS correction')
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
    cmd = 'gmt grdmath {0} {1} SUB = {2}=gd:GTiff'.format(gacos2, gacos1, outfile)
    #print(cmd)
    rc = os.system(cmd)
    if os.path.exists(outfile):
        return outfile
    else:
        print('error in GACOS processing of pair '+pair)
        return False


def remove_height_corr(ifg_ml, corr_thres = 0.5, tmpdir = os.getcwd(), dounw = True, nonlinear=False):
    '''
     first, correlate ifg_ml['pha'] and ifg_ml['hgt'] in blocks
      - better to keep dounw=True that unwraps each block by snaphu. but it can be slow
     get coefficient for correction of correlating areas
     interpolate the coefficient throughout whole raster
     multiply by hgt = 'to_remove'
    '''
    ifg_mlc = ifg_ml.copy(deep=True)
    ifg_mlc = filter_ifg_ml(ifg_mlc)
    ifg_mlc['toremove'] = 0*ifg_mlc['coh']
    t = time.process_time()
    minheight=float(ifg_ml.hgt.quantile(0.25)+200)
    thisisit, thistype = correct_hgt(ifg_mlc, blocklen = 40, tmpdir = tmpdir, dounw = dounw, nonlinear=nonlinear, minheight=minheight)
    elapsed_time = time.process_time() - t
    print('Elapsed time for hgt correction: {:f}'.format(elapsed_time))
    #thisisit can be either False, xr.DataArray, or np.float - ok, adding 'thistype' that can be bool, float, xr
    if not thistype == 'bool':
        if thistype == 'float':
            print('we use average value of {} rad/km'.format(str(thisisit*1000)))
        else:
            print('using hgt correlation grid to reduce hgt component')
        ifg_mlc['toremove'].values = thisisit*ifg_mlc['hgt']
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
            unwr = unwrap_np(cpx, coh, defomax = 0.3, tmpdir=tmpdir)
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


def unwrap_xr(ifg, mask=True, defomax = 0.3, tmpdir=os.getcwd()):
    #tmp = ifg[['cpx','coh']].copy(deep=True)
    coh = ifg.coh.values
    cpx = ifg.cpx.fillna(0).astype(np.complex64).values
    if mask:
        mask = ifg.mask.fillna(0).values
    unw = unwrap_np(cpx, coh, mask = mask, defomax = defomax, tmpdir=tmpdir, deltemp=False)
    ifg['unw']=ifg.pha.copy(deep=True)
    ifg['unw'].values = unw
    return ifg


def unwrap_np(cpx, coh, defomax = 0.3, tmpdir=os.getcwd(), mask = False, deltemp = True):
    '''
    unwraps given arrays
    '''
    if not os.path.exists(tmpdir):
        os.mkdir(tmpdir)
    try:
        binmask= os.path.join(tmpdir,'mask.bin')
        mask.astype(np.byte).tofile(binmask)
    except:
        binmask=None
    bincoh = os.path.join(tmpdir,'coh.bin')
    #binR = os.path.join(tmpdir,'R.bin')
    #binI = os.path.join(tmpdir,'I.bin')
    binCPX = os.path.join(tmpdir,'cpxifg.bin')
    unwbin = os.path.join(tmpdir,'unw.bin')
    # create R, I -> CPX as expected by snaphu
    r = np.real(cpx) #.tofile(binR)
    i = np.imag(cpx) #.tofile(binI)
    #RI2cpx(binR, binI, binCPX)
    RI2cpx(r, i, binCPX)
    #and coh
    coh.astype(np.float32).tofile(bincoh)
    # unwrap it
    width = coh.shape[0]
    #with nostdout():
    main_unwrap(binCPX, bincoh, maskbin = binmask, outunwbin = unwbin, width = width, defomax = defomax, printout = False)
    # and load it back
    unw1 = np.fromfile(unwbin,dtype=np.float32)
    unw1 = unw1.reshape(coh.shape)
    if deltemp:
        #shutil - delete tmpdir!!!!!
        shutil.rmtree(tmpdir)
    return unw1


def gausspha2cpx(ifgxr):
    ifgxr['pha'] = ifgxr.gauss_pha
    tempar_mag1 = np.ones_like(ifgxr.pha)
    cpxarr = magpha2RI_array(tempar_mag1, ifgxr.pha.values)
    ifgxr['cpx'].values = cpxarr
    return ifgxr


import dask.array as da

def correct_hgt(ifg_mlc, blocklen = 20, tmpdir = os.getcwd(), dounw = True, num_workers = 1, nonlinear=False, minheight=200, mingausscoh=0.4):
    #ifg_ml['hgtcorr'] = ifg_ml['pha']
    winsize = (blocklen, blocklen)
    cohb = da.from_array(ifg_mlc['coh'].where(ifg_mlc.gauss_coh>mingausscoh).where(ifg_mlc.hgt>minheight).astype(np.float32).fillna(0.001), chunks=winsize)
    #phab = da.from_array(ifg_ml['pha'].astype(np.float32).fillna(0), chunks=winsize)
    cpxb = da.from_array(ifg_mlc['cpx'].where(ifg_mlc.gauss_coh>mingausscoh).where(ifg_mlc.hgt>minheight).astype(np.complex64).fillna(0), chunks=winsize)
    hgtb = da.from_array(ifg_mlc['hgt'].where(ifg_mlc.gauss_coh>mingausscoh).where(ifg_mlc.hgt>minheight).astype(np.float32).fillna(0), chunks=winsize)
    f = da.map_blocks(block_hgtcorr, cpxb, cohb, hgtb, procdir = tmpdir, dounw = dounw, meta=np.array(()), chunks = (1,1))
    try:
        #with nostdout():
        hgtcorr =  f.compute(num_workers=num_workers)
    except:
        print('error in computing hgt correlation grid')
        return False, 'bool'
    # make it to xr:
    aaa = xr.Dataset()
    aaa['coh'] = ifg_mlc['coh'].fillna(0).coarsen({'lat': blocklen, 'lon': blocklen}, boundary='pad').mean()
    aaa['hgtcorr'] = aaa['coh']
    aaa.hgtcorr.values = hgtcorr
    #count means - number of non-nan data..!
    if aaa.hgtcorr.count() > 50 and nonlinear:
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
        #so this will take only median, to perform only linear heights correction
        out = np.nanmedian(hgtcorr)
        outype = 'float'
        if np.isnan(out):
            print('all NaNs in hgt corr')
            return False, 'bool'
        if np.abs(out) < 0.001:
            print('almost nothing to reduce for hgt')
            print('the estimate was: '+str(out))
            return False, 'bool'
    #if we got here, means, now it is up to splining it to the full (ml) resolution... but --- at this moment, i will just use average value
    return out, outype


def debug(frame, pair):
    #from LiCSAR_lib.unwrp_multiscale import *
    import matplotlib.pyplot as plt
    procdir = os.getcwd()
    ml = 10
    fillby = 'gauss'
    thres = 0.3
    prevest = None
    hgtcorr = True
    pre_detrend=True
    gacoscorr = True
    outtif = None
    defomax = 0.3
    add_resid = True
    smooth = True
    prev_ramp = None
    rampit=False
    cohratio = None
    keep_coh_debug = False
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
    
    
    ifg_pha_file = os.path.join(geoifgdir,pair+'.geo.diff_pha.tif')
    coh_file = os.path.join(geoifgdir,pair+'.geo.cc.tif')
    landmask_file = os.path.join(geoframedir,'metadata',frame+'.geo.landmask.tif')
    hgtfile = os.path.join(geoframedir,'metadata',frame+'.geo.hgt.tif')
    if gacoscorr:
        gacoscorrfile = os.path.join(tmpgendir,'gacos.tif')
        try:
            gacoscorrfile = make_gacos_ifg(frame, pair, gacoscorrfile)
        except:
            print('error processing gacos data for pair '+pair)
            gacoscorrfile = False
    else:
        gacoscorrfile = False
    
    
    if gacoscorrfile:
        print('GACOS data found, using to improve unwrapping')
        #ingacos = xr.open_dataset(gacoscorrfile)
        ingacos = load_tif2xr(gacoscorrfile)
    else:
        gacoscorr = False
    
    
    inpha = load_tif2xr(ifg_pha_file)
    incoh = load_tif2xr(coh_file)
    incoh = incoh/255
    inmask = incoh.copy(deep=True)
    inmask.values = np.byte(incoh > 0)
    if os.path.exists(landmask_file):
        landmask = load_tif2xr(landmask_file)
        #landmask = xr.open_dataset(landmasknc)
        inmask.values = landmask.values * inmask.values
    else:
        landmask = None
    
    
    ifg = xr.Dataset()
    ifg['pha'] = inpha
    ifg['coh'] = ifg['pha']
    ifg['coh'].values = incoh.values
    ifg['mask'] = ifg['pha']
    ifg['mask'].values = inmask.values
    ifg['mask_extent'] = ifg['pha'].where(ifg['pha'] == 0).fillna(1)
    if hgtcorr:
        if os.path.exists(hgtfile):
            try:
                inhgt = load_tif2xr(hgtfile)
                ifg['hgt'] = ifg['pha']
                ifg['hgt'].values = inhgt.values
            except:
                print('ERROR in importing heights!')
                hgtcorr = False
    
    
    if gacoscorr:
        ifg['gacos'] = ifg.pha
        ifg['gacos'].values = ingacos.values
    
    
    # erase from memory
    inpha = ''
    incoh = ''
    inmask = ''
    ingacos = ''
    inhgt = ''
    # now doing multilooking, using coh as mag...
    #make complex from coh and pha
    ifg['cpx'] = ifg.coh.copy()
    if type(cohratio) != type(None):
        ifg['cpx'].values = magpha2RI_array(cohratio.values, ifg.pha.values)
    else:
        ifg['cpx'].values = magpha2RI_array(ifg.coh.values, ifg.pha.values)
    
    
    #fixing difference in xarray version... perhaps...
    if 'lat' not in ifg.coords:
        print('warning - perhaps old xarray version - trying anyway')
        ifg = ifg.rename_dims({'x':'lon','y':'lat'})
    
    resdeg = np.abs(ifg.lat[1]-ifg.lat[0])
    if cliparea_geo:
        minclipx, maxclipx, minclipy, maxclipy = cliparea_geo.split('/')
        minclipx, maxclipx, minclipy, maxclipy = float(minclipx), float(maxclipx), float(minclipy), float(maxclipy)
        # now will clip it - lat is opposite-sorted, so need to slice from max to min in y ... plus 10 pixels on all sides
        ifg = ifg.sel(lon=slice(minclipx-10*resdeg, maxclipx+10*resdeg), lat=slice(maxclipy+10*resdeg, minclipy-10*resdeg))
    
    
    print('and now switch to multilook_normalized')
    ifg_ml = multilook_normalised(ifg, ml, tmpdir = tmpdir, hgtcorr = hgtcorr, pre_detrend = pre_detrend, prev_ramp = prev_ramp, keep_coh_debug = keep_coh_debug)


def process_ifg(frame, pair, 
        procdir = os.getcwd(), 
        ml = 10, fillby = 'gauss', 
        thres = 0.35, prevest = None, 
        hgtcorr = False, pre_detrend=True,
        gacoscorr = True, outtif = None, 
        defomax = 0.3, add_resid = True, smooth = False, 
        prev_ramp = None, rampit=False, cohratio = None, 
        keep_coh_debug = True, replace_ml_pha = None,
        coh2var = False, cliparea_geo = None,
        subtract_gacos = False, dolocal = False):
    '''
    main function for unwrapping a geocoded LiCSAR interferogram.
    ml .. multilook factor in both lat and lon directions (would use pha+coh)
    coh2var - something to try... perhaps this is better for weighting?
    cliparea_geo - as lon1/lon2/lat1/lat2
    fillby - 'nearest' or 'gauss' - gauss might be better but also slower (notice the 'gapfill iterations')
    20220506 - setting smooth to False by default - some extra errors here!
    '''
    pubdir = os.environ['LiCSAR_public']
    geoframedir = os.path.join(pubdir,str(int(frame[:3])),frame)
    if dolocal:
        geoifgdir = os.path.join('GEOC',pair)
    else:
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
    # to load orig unw_file
    #unw_file = os.path.join(geoifgdir,pair+'.geo.unw.tif')
    # do gacos if exists
    if gacoscorr:
        gacoscorrfile = os.path.join(tmpgendir,'gacos.tif')
        try:
            gacoscorrfile = make_gacos_ifg(frame, pair, gacoscorrfile)
        except:
            print('error processing gacos data for pair '+pair)
            gacoscorrfile = False
    else:
        gacoscorrfile = False
    if gacoscorrfile:
        print('GACOS data found, using to improve unwrapping')
        #ingacos = xr.open_dataset(gacoscorrfile)
        ingacos = load_tif2xr(gacoscorrfile)
    else:
        gacoscorr = False
    # prepare input dataset and load
    #use orig ifg and coh
    inpha = load_tif2xr(ifg_pha_file)
    incoh = load_tif2xr(coh_file)
    incoh = incoh/255
    inmask = incoh.copy(deep=True)
    inmask.values = np.byte(incoh > 0)
    if os.path.exists(landmask_file):
        landmask = load_tif2xr(landmask_file)
        #landmask = xr.open_dataset(landmasknc)
        inmask.values = landmask.values * inmask.values
    else:
        landmask = None
    # to get all in on encapsulement (with coords..)
    ifg = xr.Dataset()
    ifg['pha'] = inpha
    ifg['coh'] = ifg['pha']
    ifg['coh'].values = incoh.values
    ifg['mask'] = ifg['pha']
    ifg['mask'].values = inmask.values
    ifg['mask_extent'] = ifg['pha'].where(ifg['pha'] == 0).fillna(1)
    #if hgtcorr:
    # including hgt anyway - would be useful later
    if os.path.exists(hgtfile):
        try:
            inhgt = load_tif2xr(hgtfile)
            ifg['hgt'] = ifg['pha']
            ifg['hgt'].values = inhgt.values
        except:
            print('ERROR in importing heights!')
            hgtcorr = False
    if gacoscorr:
        ifg['gacos'] = ifg.pha
        ifg['gacos'].values = ingacos.values
    # erase from memory
    inpha = ''
    incoh = ''
    inmask = ''
    ingacos = ''
    inhgt = ''
    # now doing multilooking, using coh as mag...
    #make complex from coh and pha
    ifg['cpx'] = ifg.coh.copy()
    if coh2var:
        if type(cohratio) != type(None):
            # use cohratio for better weights
            coh = cohratio
        else:
            coh = ifg['coh']
        # if this is better, i will change it and have it fixed
        cohratio = (2*coh**2)/(1-coh**2)
    if type(cohratio) != type(None):
        ifg['cpx'].values = magpha2RI_array(cohratio.values, ifg.pha.values)
    else:
        ifg['cpx'].values = magpha2RI_array(ifg.coh.values, ifg.pha.values)
    #fixing difference in xarray version... perhaps...
    if 'lat' not in ifg.coords:
        print('warning - perhaps old xarray version - trying anyway')
        ifg = ifg.rename_dims({'x':'lon','y':'lat'})
    #ifg_ml = multilook_by_gauss(ifg, ml)
    # result more or less same, but much more economic way in memory..
    resdeg = np.abs(ifg.lat[1]-ifg.lat[0])
    # now crop if needed:
    if cliparea_geo:
        minclipx, maxclipx, minclipy, maxclipy = cliparea_geo.split('/')
        minclipx, maxclipx, minclipy, maxclipy = float(minclipx), float(maxclipx), float(minclipy), float(maxclipy)
        # now will clip it - lat is opposite-sorted, so need to slice from max to min in y ... plus 10 pixels on all sides
        ifg = ifg.sel(lon=slice(minclipx-10*resdeg, maxclipx+10*resdeg), lat=slice(maxclipy+10*resdeg, minclipy-10*resdeg))
        # not the best here, as pixels might get slightly shifted, but perhaps not that big deal (anyway prev_ramp is 'blurred')
        if not type(prev_ramp) == type(None):
            prev_ramp = prev_ramp.sel(lon=slice(minclipx-10*resdeg, maxclipx+10*resdeg), lat=slice(maxclipy+10*resdeg, minclipy-10*resdeg))
    #WARNING - ONLY THIS FUNCTION HAS GACOS INCLUDED NOW! (and heights fix!!!)
    ifg_ml = multilook_normalised(ifg, ml, tmpdir = tmpdir, hgtcorr = hgtcorr, pre_detrend = pre_detrend, prev_ramp = prev_ramp, keep_coh_debug = keep_coh_debug, replace_ml_pha = replace_ml_pha)
    width = len(ifg_ml.lon)
    length = len(ifg_ml.lat)
    #here we should keep the (not wrapped) phase we remove due to corrections - here, heights, and gacos
    #mask it ... this worked ok here:
    #thres = 0.5
    #gmmm.... it works better if i use the 'trick on gauss of normalised ifg....'
    #thres = 0.25
    mask_gauss = (ifg_ml.gauss_coh > thres)*1
    #return (unmask) pixels that have coh > 0.25
    mask_gauss.values[mask_gauss.values == 0] = 1*(ifg_ml.coh > 0.25).values[mask_gauss.values == 0]
    ifg_ml['mask_coh'] = mask_gauss.fillna(0)
    # additionally remove islands of size over 7x7 km
    # how many pixels are in 7x7 km region?
    lenthres = 7 # km
    # resolution of orig ifg is expected 0.1 km
    #origres = 0.1
    resdeg = np.abs(ifg.lat[1]-ifg.lat[0])
    latres = 111.32 * np.cos(np.radians(ifg_ml.lat.mean()))
    origres = float(latres * resdeg) # in km
    #
    pixels = int(round(lenthres/origres/ml))
    pixelsno = pixels**2
    npa = ifg_ml['mask_coh'].where(ifg_ml['mask_coh']==1).where(ifg_ml['mask']==1).values
    ifg_ml['mask_full'] = ifg_ml['mask_coh']
    ifg_ml['mask_full'].values = remove_islands(npa, pixelsno)
    ifg_ml['mask_full'] = ifg_ml['mask_full'].fillna(0).astype(np.int8)
    #ifg_ml['mask_coh'] = ifg_ml['mask'].where(ifg_ml.gauss_coh > thres).where(ifg_ml.coh > thres).fillna(0)
    #here the origpha is just before the gapfilling filter
    ifg_ml['origpha'] = ifg_ml.pha.copy(deep=True)
    #try without modal filter..
    #ifg_ml = filter_mask_modal(ifg_ml, 'mask_coh', 'mask_coh', 8)
    print('interpolate coh-based masked areas of gauss pha')
    #tofillpha = ifg_ml.pha.fillna(0).where(ifg_ml.mask_coh.where(ifg_ml.mask == 1).fillna(1) == 1)
    if fillby == 'gauss':
        # keep smooth always on - much better...
        if smooth:
            ifg_ml['pha'] = ifg_ml.gauss_pha
        #trying astropy approach now:
        print('filling through Gaussian kernel')
        kernel = Gaussian2DKernel(x_stddev=2)
        # create a "fixed" image with NaNs replaced by interpolated values
        tempar_mag1 = np.ones_like(ifg_ml.pha)
        #cpxarr = magpha2RI_array(tempar_mag1, ifg_ml.pha.fillna(0).values)
        cpxarr = magpha2RI_array(tempar_mag1, ifg_ml.pha.values)
        ifg_ml['cpx_tofill'] = ifg['cpx']
        ifg_ml['cpx_tofill'].values = cpxarr
        # using only extent - avoiding landmask as rivers would have problems
        #tofill = ifg_ml['cpx_tofill'].where(ifg_ml.mask_coh.where(ifg_ml.mask == 1).fillna(1) == 1)
        #tofill = ifg_ml['cpx_tofill'].where(ifg_ml.mask_coh.where(ifg_ml.mask_extent == 1).fillna(1) == 1)
        # or not... let's mask fully
        i = 1
        tofill = ifg_ml['cpx_tofill'].where(ifg_ml.mask_full.where(ifg_ml.mask_extent == 1).fillna(1) == 1)
        tofillR = np.real(tofill)
        tofillI = np.imag(tofill)
        filledR = interpolate_replace_nans(tofillR.values, kernel)
        filledI = interpolate_replace_nans(tofillI.values, kernel)
        ifg_ml['pha'].values = np.angle(filledR + 1j*filledI)
        #sometimes the whole area is not within gauss kernel
        while np.max(np.isnan(ifg_ml['pha'].values)):
            i = i+1
            print('gapfilling iteration '+str(i))
            if i>3:
                # no need to add more heavy iterations
                print('filling by nearest neighbours')
                #tofillpha = ifg_ml.pha.where(ifg_ml.mask_coh.where(ifg_ml.mask == 1).fillna(1) == 1)
                #tofillpha = ifg_ml.pha.where(ifg_ml.mask_full.where(ifg_ml.mask_extent == 1).fillna(1) == 1)
                ifg_ml['pha'].values = interpolate_nans(ifg_ml['pha'].values, method='nearest')
            else:
                cpxarr = magpha2RI_array(tempar_mag1, ifg_ml.pha.values)
                ifg_ml['cpx_tofill'].values = cpxarr
                tofill = ifg_ml['cpx_tofill']
                tofillR = np.real(tofill)
                tofillI = np.imag(tofill)
                filledR = interpolate_replace_nans(tofillR.values, kernel)
                filledI = interpolate_replace_nans(tofillI.values, kernel)
                ifg_ml['pha'].values = np.angle(filledR + 1j*filledI)
        #sometimes the whole area is not within gauss kernel - use NN for that:
        #if np.max(np.isnan(ifg_ml['pha'].values)):
        #    # no, this would be too far. using only filling by zero
        #    ifg_ml['pha'] = ifg_ml['pha'].fillna(0)
        #    #ifg_ml['pha'].values = interpolate_nans(ifg_ml['pha'].values, method='nearest')
        #ifg_ml.pha.fillna(0).where(ifg_ml.mask_coh.where(ifg_ml.mask == 1).fillna(1) == 1)
        #ifg_ml['pha'].values = interpolate_replace_nans(tofillpha.values, kernel)
    elif fillby == 'none':
        print('skipping any nan filling')
        ifg_ml['pha'] = ifg_ml.gauss_pha
    else:
        print('filling by nearest neighbours')
        #tofillpha = ifg_ml.pha.fillna(0).where(ifg_ml.mask_coh.where(ifg_ml.mask == 1).fillna(1) == 1)
        #tofillpha = ifg_ml.pha.fillna(0).where(ifg_ml.mask_coh.where(ifg_ml.mask == 1).fillna(1) == 1)
        tofillpha = ifg_ml.pha.where(ifg_ml.mask_full.where(ifg_ml.mask_extent == 1).fillna(1) == 1)
        ifg_ml['pha'].values = interpolate_nans(tofillpha.values, method='nearest')
        #ifg_ml['gauss_pha'] = ifg_ml['gauss_pha'].fillna(0)
    print('debug: now pha is fine-filled layer but with some noise at edges - why is that? not resolved. so adding one extra gauss filter')
    # ok, i see some high freq signal is still there.. so filtering once more (should also help after the nan filling)
    if smooth:
        print('an extra Gaussian smoothing here')
        ifg_ml = filter_ifg_ml(ifg_ml)
        ifg_ml['pha'] = ifg_ml['gauss_pha']
    #exporting for snaphu
    #normalise mag from the final pha
    tempar_mag1 = np.ones_like(ifg_ml.pha)
    cpxarr = magpha2RI_array(tempar_mag1, ifg_ml.pha.fillna(0).values) #no need to fillna, but just in case...
    ifg_ml['gauss_cpx'].values = cpxarr
    print('unwrapping by snaphu')
    binmask= os.path.join(tmpdir,'gaussmask.bin')
    #bincoh = os.path.join(tmpdir,'gausscoh.bin')
    bincoh = os.path.join(tmpdir,'coh.bin')
    #binR = os.path.join(tmpdir,'gaussR.bin')
    #binI = os.path.join(tmpdir,'gaussI.bin')
    binCPX = os.path.join(tmpdir,'cpxgaussifg.bin')
    outunwbin = os.path.join(tmpdir,'gaussunwrapped.bin')
    #print('exporting to bin files')
    #ifg_ml.mask_coh.fillna(0).values.astype(np.byte).tofile(binmask)
    # full masking may be too much for snaphu here:
    #ifg_ml.mask_full.fillna(0).values.astype(np.byte).tofile(binmask)
    ifg_ml.mask.fillna(0).values.astype(np.byte).tofile(binmask)
    #ifg_ml.gauss_coh.fillna(0.001).values.astype(np.float32).tofile(bincoh)
    # we should use the orig coh for weights... and perhaps very low coh values instead of 0
    ifg_ml.coh.fillna(0.001).values.astype(np.float32).tofile(bincoh)
    r = np.real(ifg_ml.gauss_cpx).astype(np.float32).fillna(0).values #.tofile(binR)
    i = np.imag(ifg_ml.gauss_cpx).astype(np.float32).fillna(0).values #.tofile(binI)
    RI2cpx(r, i, binCPX)
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
        #print('unwrapping')
        #main_unwrap(binCPX, bincoh, binmask, outunwbin, width, defomax = defomax, printout=False)
        # 2022-01-14 - avoiding mask here - it does only worse
        #main_unwrap(binCPX, bincoh, None, outunwbin, width, defomax = defomax, printout=False)
        # 2022-04-04 - returning the mask! result is really bad with it, at least at islands!
        main_unwrap(binCPX, bincoh, binmask, outunwbin, width, defomax = defomax, printout=False)
    print('importing snaphu result to ifg_ml')
    binfile = outunwbin
    #toxr = ifg_ml
    daname = 'unw'
    dtype = np.float32
    unw1 = np.fromfile(binfile,dtype=dtype)
    unw1 = unw1.reshape(ifg_ml.pha.shape)
    #unw1 = np.flip(unw1,axis=0)
    ifg_ml[daname] = ifg_ml['pha'] #.copy(deep=True)
    ifg_ml[daname].values = unw1
    #ok, so the gauss-based coh mask is not the best to do... so exporting 'all pixels'
    #ifg_ml[daname] = ifg_ml[daname]*ifg_ml['mask']
    #print('20210722 - testing now - using gauss-based coh mask, ignore the next message:')
    #ifg_ml[daname] = ifg_ml[daname]*ifg_ml['mask_coh']
    #this would do median correction, not doing it now:
    #ifg_ml[daname].values = ifg_ml[daname].values - np.nanmedian(ifg_ml[daname].where(ifg_ml.mask_coh>0).values)
    #ifg_ml[daname] = ifg_ml[daname]*ifg_ml.mask_coh
    #ifg_ml[daname] = ifg_ml[daname]*ifg_ml.mask
    #ifg_ml[daname] = ifg_ml[daname]*ifg_ml['mask_coh']
    ifg_ml[daname] = ifg_ml[daname]*ifg_ml['mask_full']
    #print('unwrap also residuals from the filtered cpx, and add to the final unw - mask only waters..')
    cpxarr = magpha2RI_array(tempar_mag1, ifg_ml.origpha.fillna(0).where(ifg_ml.mask == 1).values)
    ifg_ml['origcpx'] = ifg_ml['gauss_cpx'] #.copy(deep=True)
    ifg_ml['origcpx'].values = cpxarr
    if add_resid:
        print('unwrapping residuals and adding back to the final unw output')
        # maybe can use this somehow? by yma. i checked and it is correct
        # delta = np.angle(np.exp(np.complex(0+1j)*( np.angle(ifg_filt) - np.angle(ifg_unfilt))))
        # delta = np.angle(np.exp(0+1j)*( np.angle(ifg_filt) - np.angle(ifg_unfilt)))
        # delta = np.angle(np.exp(0+1j)*( pha_filt - pha_unfilt ) )
        # ifg_unw_unfilt = ifg_unw_filt - delta
        ifg_ml['resid_cpx'] = ifg_ml.origcpx * ifg_ml.gauss_cpx.conjugate() #* ifg_ml.mask_full
        incpx = 'resid_cpx'
        #binR = os.path.join(tmpdir,incpx+'.R.bin')
        #binI = os.path.join(tmpdir,incpx+'.I.bin')
        binCPX = os.path.join(tmpdir,incpx+'.cpx.bin')
        outunwbin = os.path.join(tmpdir,incpx+'.unw.bin')
        binfile = outunwbin
        daname = incpx+'.unw'
        #
        r = np.real(ifg_ml[incpx]).astype(np.float32).fillna(0).values #.tofile(binR)
        i = np.imag(ifg_ml[incpx]).astype(np.float32).fillna(0).values #.tofile(binI)
        RI2cpx(r, i, binCPX)
        #main_unwrap(binCPX, bincoh, binmask, outunwbin, width, defomax = defomax/2)
        # ok, just hold the defomax down - discontinuities are not wanted or expected here
        main_unwrap(binCPX, bincoh, binmask, outunwbin, width, defomax = 0, printout = False)
        unw1 = np.fromfile(binfile,dtype=dtype)
        unw1 = unw1.reshape(ifg_ml.pha.shape)
        ifg_ml[daname] = ifg_ml['pha'] #.copy()
        ifg_ml[daname].values = unw1
        #ifg_ml[daname] = ifg_ml[daname]*ifg_ml['mask']
        #ifg_ml[daname] = ifg_ml[daname]*ifg_ml['mask_coh']
        # ensure values are correctly surrounded by zeroes - i.e. shifting by points masked away
        # 2022/02 - actually seems weird to me. skipping. also, i added mask binary to update snaphu cost solution..
        #ifg_ml[daname] = ifg_ml[daname] - np.nanmedian(ifg_ml[daname].where(ifg_ml.mask_extent==0).values)
        #ifg_ml[daname] = ifg_ml[daname]*ifg_ml['mask_full']
        #ifg_ml[daname] = ifg_ml[daname] - np.nanmedian(ifg_ml[daname].where(ifg_ml.mask_extent==1).values)
        #ifg_ml[daname] = ifg_ml[daname]*ifg_ml['mask_full']
        #nanmed = np.nanmedian(ifg_ml[daname].where(ifg_ml.mask_coh>0).values)
        nanmed = np.nanmedian(ifg_ml[daname].where(ifg_ml.mask_full>0).values)
        ifg_ml[daname].values = ifg_ml[daname].values - nanmed
        ifg_ml[daname] = ifg_ml[daname]*ifg_ml['mask_full']
        #print('debug - avoiding median correction now, although maybe ok for add_resid: median was {} rad'.format(str(nanmed)))
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
    #ifg_ml['unw'] = ifg_ml['unw'].where(ifg_ml.mask>0)
    # ok ok.... let's mask by the gauss mask.. although, can be quite missing lot of areas
    #ifg_ml['unw'] = ifg_ml['unw'].where(ifg_ml.mask_coh>0)
    # hmmm... just to make it nicer...
    #ifg_ml['unw'] = ifg_ml['unw'].where(ifg_ml.mask>0)
    ifg_ml['unw'] = ifg_ml['unw'].where(ifg_ml.mask_full>0)
    # but this may be wrong!
    # 2022-04-04 - ok, but gacos and some other corrections based on model would bring full shift,
    # including the one already inside range offsets during coreg. so we should indeed shift by median:
    ifg_ml['unw'] = ifg_ml['unw'] - ifg_ml['unw'].median()
    if rampit:
        ifg_ml['unw_orig'] = ifg_ml['unw'].copy()
        ifg_ml['unw'].values= interpolate_replace_nans(ifg_ml['unw'].values, kernel)
        ifg_ml['unw'].values = filter_nan_gaussian_conserving(ifg_ml['unw'].values, sigma=2, trunc=4)
        ifg_ml['unw'] = ifg_ml['unw'].fillna(0).where(ifg_ml.mask_extent>0)
    #else:
    #    ifg_ml['unw'] = ifg_ml['unw'].where(ifg_ml.mask_full>0)
    # finally clip again, without border pixels:
    if cliparea_geo:
        ifg_ml = ifg_ml.sel(lon=slice(minclipx, maxclipx), lat=slice(maxclipy, minclipy))
    if add_resid:
        # 2022-04-04 - fixing for final residuals - without unwrapping them - in case the unw had some spatially propagating error
        # BUT ALSO the 'weird vertical lines' sometimes induced after snaphu unwrapping. those lines have, however, values very close to 0, but not 0
        # it generally should not be needed though!
        if 'origpha_noremovals' in ifg_ml:
            residpha = wrap2phase(wrap2phase(ifg_ml['unw'].fillna(0).values) - ifg_ml.origpha_noremovals)
        else:
            tempar_mag1 = np.ones_like(ifg_ml.pha)
            cpxarr = magpha2RI_array(tempar_mag1, wrap2phase(ifg_ml['unw'].fillna(0).values))
            ifg_ml['unwcpx'] = ifg_ml['gauss_cpx']
            ifg_ml['unwcpx'].values = cpxarr
            ifg_ml['resid_final'] = ifg_ml.unwcpx * ifg_ml.origcpx.conjugate() #* ifg_ml.mask_full
            # but ignoring the overall shift as we do the median shifting before
            residpha = np.angle(ifg_ml['resid_final'].where(ifg_ml.mask_full>0))
        medres = np.nanmedian(residpha)
        residpha = residpha - medres
        # ifg_ml['resid_final'].values = residpha; ifg_ml['resid_final'].plot(); plt.show()
        print('final check for residuals: their std is '+str(np.nanstd(residpha)))
        # ok, so i assume that the unw would not help anymore, so just adding to unw as it is
        ifg_ml['unw'] = ifg_ml['unw'] - residpha
    # now, we may need to save without gacos itself:
    if subtract_gacos:
        if 'gacos' in ifg_ml:
            ifg_ml['unw'] = ifg_ml['unw'] - (ifg_ml['gacos'] - ifg_ml['gacos'].median())
            ifg_ml['unw'] = ifg_ml['unw'].where(ifg_ml.mask_full>0)
    if outtif:
        #ifg_ml.pha.fillna(0).where(ifg_ml.mask_coh.where(ifg_ml.mask == 1).fillna(1) == 1)
        #ifg_ml['unw'].to_netcdf(outtif+'.nc')
        toexp = ifg_ml['unw'].copy(deep=True)
        export_xr2tif(toexp, outtif)
        #toexp.values = np.flipud(toexp.values)
        # this does not work, probably due to multilooking... yikes
        #ifg_ml['unw'].to_netcdf(outtif+'.nc')
        #rc = os.system('gmt grdconvert -G{0}=gd:GTiff -R{1} {0}.nc'.format(outtif, ifg_pha_file))
        #rc = os.system('source {0}/lib/LiCSAR_bash_lib.sh; create_preview_unwrapped {1} {2}'.format(os.environ['LiCSARpath'], outtif, frame))
        rc = os.system('source {0}/lib/LiCSAR_bash_lib.sh; create_preview_unwrapped {1}'.format(os.environ['LiCSARpath'], outtif))
        #try:
        #    os.remove(outtif+'.nc')
        #except:
        #    print('ERROR removing the nc file - something wrong with export')
    return ifg_ml



def export_xr2tif(xrda, tif, lonlat = True, debug = True, dogdal = True):
    import rioxarray
    #coordsys = xrda.crs.split('=')[1]
    if debug:
        xrda = xrda.astype(np.float32)
    coordsys = "epsg:4326"
    if lonlat:
        xrda = xrda.rio.set_spatial_dims(x_dim="lon", y_dim="lat", inplace=True)
    else:
        xrda = xrda.rio.set_spatial_dims(x_dim="x", y_dim="y", inplace=True)
    xrda = xrda.rio.write_crs(coordsys, inplace=True)
    if dogdal:
        xrda.rio.to_raster(tif+'tmp.tif', compress='deflate')
        cmd = 'gdalwarp -t_srs EPSG:4326 {0} {1}'.format(tif+'tmp.tif', tif)
        runcmd(cmd, printcmd = False)
        os.remove(tif+'tmp.tif')
    else:
        xrda.rio.to_raster(tif, compress='deflate')


def process_frame(frame, ml = 10, hgtcorr = True, cascade=False, use_amp_stab = False,
            use_coh_stab = False, keep_coh_debug = True, export_to_tif = False, 
            gacoscorr = True, phase_bias_experiment = False, cliparea_geo = None,
            pairsetfile = None, subtract_gacos = False, nproc = 1, smooth = False, dolocal = False):
    '''
    hint - try use_coh_stab = True.. maybe helps against loop closure errors?!
    '''
    #if cascade and ml>1:
    #    print('error - the cascade approach is ready only for ML1')
    #    return False
    #the best to run in directory named by the frame id
    pubdir = os.environ['LiCSAR_public']
    geoframedir = os.path.join(pubdir,str(int(frame[:3])),frame)
    if dolocal:
        geoifgdir = 'GEOC'
        if not os.path.exists(geoifgdir):
            print('ERROR: the GEOC directory does not exist, cancelling')
            exit()
    else:
        geoifgdir = os.path.join(geoframedir,'interferograms')
    inputifgdir = geoifgdir
    if phase_bias_experiment:
        print('running for the phas bias experiment - make sure you are inside your folder with ifgs, e.g.')
        print('/work/scratch-pw/earyma/LiCSBAS/bias/138D_05142_131313/wrapped/GEOC_wrapped_ml10')
        inputifgdir = os.getcwd()
    hgtfile = os.path.join(geoframedir,'metadata', frame+'.geo.hgt.tif')
    raster = gdal.Open(hgtfile)
    framewid = raster.RasterXSize
    framelen = raster.RasterYSize
    #if cliparea_geo:
    #    import rioxarray as rio
    #    hgt = xr.open_dataarray(hgtfile)
    cohratio = None
    if use_amp_stab:
        print('calculating amplitude stability')
        try:
            ampstabfile = frame+'_ampstab.nc'
            if os.path.exists(ampstabfile):
                print('using existing ampstabfile')
                ampstab = xr.open_dataarray(ampstabfile)
            else:
                ampavg, ampstd = build_amp_avg_std(frame)
                ampstab = 1 - ampstd/ampavg
                ampstab.values[ampstab<=0] = 0.00001
                ampstab.to_netcdf(ampstabfile)
        except:
            print('some error happened, disabling use of amplitude stability')
            use_amp_stab = False
    if use_coh_stab:
        cohstabdays = 12
        print('calculating coherence stability, using only {} days coherences'.format(str(cohstabdays)))
        try:
            cohratiofile = frame+'_cohratio.nc'
            if os.path.exists(cohratiofile):
                print('using existing cohratiofile')
                cohratio = xr.open_dataarray(cohratiofile)
            else:
                cohavg, cohstd = build_coh_avg_std(frame, days = cohstabdays, monthly = False)
                # ok, original coh_stab = 1 - coh_dispersion, i.e.:
                cohratio = 1 - cohstd/cohavg
                cohratio.values[cohratio<=0] = 0.00001
                print('storing DQ=1-cohstd/cohavg to '+cohratiofile)
                cohratio.to_netcdf(cohratiofile)
                # but i want now to have it logarithmic, so:
                #cohratio = cohavg/cohstd
                #hmm... not really any difference. so using the orig way
        except:
            print('some error happened, disabling use of coh stability')
            use_coh_stab = False
    pairset = None
    if pairsetfile:
        try:
            pairset = pd.read_csv(pairsetfile, header = None)
            pairset = list(pairset[0].values)
        except:
            print('error loading pairset, doing all')
            pairset = None
    if not pairset:
        pairset = os.listdir(inputifgdir)
    '''
    for multiprocessing - i need to make the below code as a function (but check, as it must use global variables! i need to learn with them better)
    then the situation would be super easy, i.e.:
    def check_and_process_ifg(pair):
        ...
    if nproc>1:
        try:
            from pathos.multiprocessing import ProcessingPool as Pool
        except:
            print('pathos not installed - not parallelism')
            nproc = 1
    if nproc>1:
        p = Pool(nproc)
        outs = p.map(check_and_process_ifg, pairset)  # out is one output per pair -> list
        p.close()  # or not?
    else:
        for pair in pairset:
            out=check_and_process_ifg(pair)
    '''
    #so try this:
    def check_and_process_ifg(pair):
        if not os.path.exists(os.path.join(geoifgdir, pair, pair+'.geo.diff_pha.tif')):
            return False
        else:
            #check its dimensions..
            raster = gdal.Open(os.path.join(geoifgdir, pair, pair+'.geo.diff_pha.tif'))
            if (framewid != raster.RasterXSize) or (framelen != raster.RasterYSize):
                #use tolerance of max pixels
                maxpixels = 4
                if ((abs(framewid - raster.RasterXSize) > maxpixels) or (abs(framelen - raster.RasterYSize) > maxpixels)):
                    print('ERROR - the file {} has unexpected dimensions, skipping'.format(os.path.join(geoifgdir, pair, pair+'.geo.diff_pha.tif')))
                    return False
                else:
                    print('ERROR - the file {} has unexpected dimensions, trying to fix'.format(os.path.join(geoifgdir, pair, pair+'.geo.diff_pha.tif')))
                    for tif in glob.glob(os.path.join(geoifgdir, pair, pair+'.geo.*.tif')):
                        outfile = tif+'.tmp.tif'
                        try:
                            filedone = reproject_to_match(tif, hgtfile, outfile)
                            if os.path.exists(outfile):
                                shutil.move(outfile, tif)
                        except:
                            print('something wrong during reprojection, skipping')
                            #continue
                            return False
                    #os.system('gmt grdsample {0} -G{1}')
            try:
                raster = gdal.Open(os.path.join(geoifgdir, pair, pair+'.geo.diff_pha.tif'))
                if (framewid != raster.RasterXSize) or (framelen != raster.RasterYSize):
                    print('ERROR - the file {} has unexpected dimensions, skipping'.format(os.path.join(geoifgdir, pair, pair+'.geo.diff_pha.tif')))
                    #continue
                    return False
            except:
                print('some error processing file {}'.format(os.path.join(geoifgdir, pair, pair+'.geo.diff_pha.tif')))
                #continue
                return False
            if not os.path.exists(os.path.join(pair,pair+'.unw')):
                print('processing pair '+pair)
                if export_to_tif:
                    outtif = os.path.join(pair,pair+'.geo.unw.tif')
                else:
                    outtif = None
                try:
                    if cascade:
                        ifg_ml = cascade_unwrap(frame, pair, downtoml = ml, procdir = os.getcwd(), outtif = outtif, subtract_gacos = subtract_gacos, smooth = smooth, hgtcorr = hgtcorr, cliparea_geo = cliparea_geo, dolocal=dolocal)
                    else:
                        #ifg_ml = process_ifg(frame, pair, procdir = os.getcwd(), ml = ml, hgtcorr = hgtcorr, fillby = 'gauss')
                        defomax = 0.3
                        replace_ml_pha = None
                        if phase_bias_experiment:
                            replace_ml_pha = os.path.join(pair, pair+'.diff_pha_cor')
                        ifg_ml = process_ifg(frame, pair, procdir = os.getcwd(), ml = ml, hgtcorr = hgtcorr, fillby = 'gauss', 
                                 thres = 0.3, defomax = defomax, add_resid = True, outtif = outtif, cohratio = cohratio, smooth = smooth,
                                 keep_coh_debug = keep_coh_debug, gacoscorr = gacoscorr, replace_ml_pha = replace_ml_pha, cliparea_geo = cliparea_geo,
                                 subtract_gacos = subtract_gacos, dolocal = dolocal)
                    (ifg_ml.unw.where(ifg_ml.mask_full > 0).values).astype(np.float32).tofile(pair+'/'+pair+'.unw')
                    ((ifg_ml.coh.where(ifg_ml.mask > 0)*255).astype(np.byte).fillna(0).values).tofile(pair+'/'+pair+'.cc')
                    # export 
                    width = len(ifg_ml.lon)
                    create_preview_bin(pair+'/'+pair+'.unw', width, ftype = 'unw')
                    os.system('rm '+pair+'/'+pair+'.unw.ras')
                    os.system('rm -r '+pair+'/'+'temp_'+str(ml))
                    os.system('rm -r '+pair+'/'+'temp_gen')
                except:
                    print('ERROR processing of pair '+pair)
                    os.system('rm -r '+pair)
            if not os.path.exists(os.path.join(pair,pair+'.unw')):
                print('some error occured and the unw was not processed')
                os.system('rm -r '+pair)
                return False
            else:
                return True
    def fix_additionals():
        mlipar = 'slc.mli.par'
        if not os.path.exists(mlipar):
            f = open(mlipar, 'w')
            f.write('range_samples: '+str(len(ifg_ml.lon))+'\n')
            f.write('azimuth_lines: '+str(len(ifg_ml.lat))+'\n')
            f.write('radar_frequency: 5405000000.0 Hz\n')
            f.close()
        if not os.path.exists('hgt'):
            if 'hgt' in ifg_ml:
                ifg_ml['hgt'].astype(np.float32).values.tofile('hgt')
        if not os.path.exists('EQA.dem_par'):
            post_lon=np.round(float(ifg_ml.lon[1] - ifg_ml.lon[0]),6)
            post_lat=np.round(float(ifg_ml.lat[1] - ifg_ml.lat[0]),6)
            cor_lat = np.round(float(ifg_ml.lat[0]),6)
            cor_lon = np.round(float(ifg_ml.lon[0]),6)
            create_eqa_file('EQA.dem_par',len(ifg_ml.lon),len(ifg_ml.lat),cor_lat,cor_lon,post_lat,post_lon)
    if nproc>1:
        try:
            from pathos.multiprocessing import ProcessingPool as Pool
        except:
            print('pathos not installed - not parallelism')
            nproc = 1
    if nproc>1:
        try:
            p = Pool(nproc)
            outs = p.map(check_and_process_ifg, pairset)  # out is one output per pair -> list
            p.close()  # or not?
            fix_additionals()
        except:
            print('some error appeared - please try manually (debug). now, just returning to no parallelism')
            nproc = 1
    if nproc == 1:
        for pair in pairset:
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
                    if export_to_tif:
                        outtif = os.path.join(pair,pair+'.geo.unw.tif')
                    else:
                        outtif = None
                    try:
                        if cascade:
                            #ifg_ml = cascade_unwrap(frame, pair, procdir = os.getcwd(), outtif = outtif, subtract_gacos = subtract_gacos, cliparea_geo = cliparea_geo)
                            ifg_ml = cascade_unwrap(frame, pair, downtoml = ml, procdir = os.getcwd(), outtif = outtif, subtract_gacos = subtract_gacos, smooth = smooth, cliparea_geo = cliparea_geo, dolocal = dolocal)
                        else:
                            #ifg_ml = process_ifg(frame, pair, procdir = os.getcwd(), ml = ml, hgtcorr = hgtcorr, fillby = 'gauss')
                            defomax = 0.3
                            replace_ml_pha = None
                            if phase_bias_experiment:
                                replace_ml_pha = os.path.join(pair, pair+'.diff_pha_cor')
                            ifg_ml = process_ifg(frame, pair, procdir = os.getcwd(), ml = ml, hgtcorr = hgtcorr, fillby = 'gauss', 
                                     thres = 0.3, defomax = defomax, add_resid = True, outtif = outtif, cohratio = cohratio, 
                                     keep_coh_debug = keep_coh_debug, gacoscorr = gacoscorr, replace_ml_pha = replace_ml_pha, cliparea_geo = cliparea_geo,
                                     subtract_gacos = subtract_gacos, dolocal = dolocal, smooth = smooth)
                        #else:
                        
                        #    print('ML set to 1 == will try running the cascade (multiscale) unwrapping')
                        #    ifg_ml = cascade_unwrap(frame, pair, procdir = os.getcwd())
                        #ifg_ml.unw.values.tofile(pair+'/'+pair+'.unw')
                        #np.flipud(ifg_ml.unw.fillna(0).values).tofile(pair+'/'+pair+'.unw')
                        #np.flipud(ifg_ml.unw.where(ifg_ml.mask_coh > 0).values).astype(np.float32).tofile(pair+'/'+pair+'.unw')
                        #np.flipud(ifg_ml.unw.where(ifg_ml.mask > 0).values).astype(np.float32).tofile(pair+'/'+pair+'.unw')
                        #hey, it works ok without flipping - probably thanks to new sort of lat values when loading from tiffs directly
                        #np.flipud(ifg_ml.unw.where(ifg_ml.mask_full > 0).values).astype(np.float32).tofile(pair+'/'+pair+'.unw')
                        (ifg_ml.unw.where(ifg_ml.mask_full > 0).values).astype(np.float32).tofile(pair+'/'+pair+'.unw')
                        #or use gauss coh here?
                        #np.flipud((ifg_ml.coh.where(ifg_ml.mask > 0)*255).astype(np.byte).fillna(0).values).tofile(pair+'/'+pair+'.cc')
                        ((ifg_ml.coh.where(ifg_ml.mask > 0)*255).astype(np.byte).fillna(0).values).tofile(pair+'/'+pair+'.cc')
                        
                        # export 
                        mlipar = 'slc.mli.par'
                        if not os.path.exists(mlipar):
                            f = open(mlipar, 'w')
                            f.write('range_samples: '+str(len(ifg_ml.lon))+'\n')
                            f.write('azimuth_lines: '+str(len(ifg_ml.lat))+'\n')
                            f.write('radar_frequency: 5405000000.0 Hz\n')
                            f.close()
                        if not os.path.exists('hgt'):
                            if 'hgt' in ifg_ml:
                                ifg_ml['hgt'].astype(np.float32).values.tofile('hgt')
                        if not os.path.exists('EQA.dem_par'):
                            post_lon=np.round(float(ifg_ml.lon[1] - ifg_ml.lon[0]),6)
                            post_lat=np.round(float(ifg_ml.lat[1] - ifg_ml.lat[0]),6)
                            cor_lat = np.round(float(ifg_ml.lat[0]),6)
                            cor_lon = np.round(float(ifg_ml.lon[0]),6)
                            create_eqa_file('EQA.dem_par',len(ifg_ml.lon),len(ifg_ml.lat),cor_lat,cor_lon,post_lat,post_lon)
                        width = len(ifg_ml.lon)
                        create_preview_bin(pair+'/'+pair+'.unw', width, ftype = 'unw')
                        #os.system('rm -r '+pair+'/'+'temp_'+str(ml))
                        os.system('rm -r '+pair+'/'+'temp_*')
                    except:
                        print('ERROR processing of pair '+pair)
                        os.system('rm -r '+pair)
                if not os.path.exists(os.path.join(pair,pair+'.unw')):
                    print('some error occured and the unw was not processed')
                    os.system('rm -r '+pair)


def create_eqa_file(eqafile,wid,nlines,cor_lat,cor_lon,post_lat,post_lon):
    f = open(eqafile, 'w')
    f.write('data_format:        REAL*4\n')
    f.write('DEM_hgt_offset:          0.00000\n')
    f.write('DEM_scale:               1.00000\n')
    f.write('width: '+str(wid)+'\n')
    f.write('nlines: '+str(nlines)+'\n')
    f.write('corner_lat: '+str(cor_lat)+'  decimal degrees\n')
    f.write('corner_lon: '+str(cor_lon)+'  decimal degrees\n')
    f.write('post_lat: '+str(post_lat)+' decimal degrees\n')
    f.write('post_lon: '+str(post_lon)+' decimal degrees\n')
    f.write('ellipsoid_name: WGS 84\n')
    f.write('ellipsoid_ra:        6378137.000   m\n')
    f.write('ellipsoid_reciprocal_flattening:  298.2572236\n')
    f.write('datum_name: WGS 1984\n')
    f.close()


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


def get_fft_std(inarr):
    '''
    a is numpy array
    improvised way, not reading much about it... pure intuition (knowing this is only first step to do it right)
    '''
    a=inarr.copy()
    a[np.isnan(a)]=0
    f = np.fft.fft2(a)
    fshift = np.fft.fftshift(f)
    magnitude_spectrum = 20*np.log(np.abs(fshift))
    return magnitude_spectrum.mean()


def multilook_normalised(ifg, ml = 10, tmpdir = os.getcwd(), hgtcorr = True, pre_detrend = True, prev_ramp = None, thres_pxcount = None, keep_coh_debug = True, replace_ml_pha = None):
    '''
    prev_ramp is an xarray dataframe, it can be of different multilooking
    '''
    #landmask it and multilook it
    if ml > 1:
        # that's for multilooking - in case of cohratio, we want to only weight phases based on that, and then return to coh
        bagcpx = ifg[['cpx']].where(ifg.mask>0).coarsen({'lat': ml, 'lon': ml}, boundary='trim')
        ifg_ml = bagcpx.sum() / bagcpx.count()   # this really equals np.nanmean, or, bagcpx.mean() - same result
        # if we use coh instead of amplitude, it may get underestimated after ML (as found by Jack McG.), so just averaging it here:
        #if type(ml_weights) != type(None):
        bagcoh = ifg[['coh']].where(ifg.mask>0).coarsen({'lat': ml, 'lon': ml}, boundary='trim')
        #coh_ml = bagcoh.sum() / bagcoh.count()
        coh_ml = bagcoh.mean()
        ifg_ml['coh'] = ifg_ml.cpx
        ifg_ml['coh'].values = coh_ml.coh.values
        # non-nan px count per window
        ifg_ml['pxcount'] = ifg_ml.cpx
        ifg_ml['pxcount'].values = bagcpx.count().cpx.values.astype(np.uint16)    # to use later - evaluate bad ML data, e.g. mask those pxls
    else:
        ifg_ml = ifg[['cpx']].where(ifg.mask>0)
        ifg_ml['coh'] = ifg_ml.cpx
        ifg_ml['coh'].values = ifg.coh.values
    #downsample mask
    if ml > 1:
        ifg_ml['mask'] = ifg.mask.coarsen({'lat': ml, 'lon': ml}, boundary='trim').max().astype(np.int8)
        ifg_ml['mask_extent'] = ifg.mask_extent.coarsen({'lat': ml, 'lon': ml}, boundary='trim').max().astype(np.int8)
        if 'gacos' in ifg.variables:
            ifg_ml['gacos'] = ifg.gacos.coarsen({'lat': ml, 'lon': ml}, boundary='trim').mean()  # or median?
        if 'hgt' in ifg.variables:
            ifg_ml['hgt'] = ifg.hgt.coarsen({'lat': ml, 'lon': ml}, boundary='trim').mean()
    else:
        ifg_ml['mask'] = ifg.mask
        ifg_ml['mask_extent'] = ifg.mask_extent
        if 'gacos' in ifg.variables:
            ifg_ml['gacos'] = ifg.gacos
        if 'hgt' in ifg.variables:
            ifg_ml['hgt'] = ifg.hgt
    #keep the original original pha values
    ifg_ml['origpha_noremovals'] = ifg_ml.cpx #.copy(deep=True)
    ifg_ml['origpha_noremovals'].values = np.angle(ifg_ml.cpx) #.copy(deep=True)
    #prepare 'toremove' layer
    ifg_ml['toremove'] = ifg_ml.cpx
    ifg_ml['toremove'].values = 0*np.angle(ifg_ml.cpx) # just make them zeroes
    if replace_ml_pha:
        # a quick fix to load other existing phase - e.g. after bias correction...
        ifg_ml['origpha_noremovals'].values = np.fromfile(replace_ml_pha, dtype=np.float32).reshape(ifg_ml['origpha_noremovals'].values.shape)
        
        cpxa = magpha2RI_array(ifg_ml.coh.values, ifg_ml.origpha_noremovals.values)
        ifg_ml['cpx'].values = cpxa
    if keep_coh_debug:
        #ok, return coh, phase back to cpx
        cpxa = magpha2RI_array(ifg_ml.coh.values, ifg_ml.origpha_noremovals.values)
        ifg_ml['cpx'].values = cpxa
    else:
        #print('debug: trying to use the cohratio rather than current coh. maybe wrong?')
        ifg_ml['orig_coh'] = ifg_ml.coh.copy(deep=True)
        ifg_ml.coh.values = np.abs(ifg_ml['cpx'])
    #remove previous estimates
    if not type(prev_ramp) == type(None):
        width = len(ifg_ml.lon)
        length = len(ifg_ml.lat)
        #prev_width = len(prev_ramp.lon)
        #prev_length = len(prev_ramp.lat)
        prev_ramp = prev_ramp.fillna(0)
        resized = cv2.resize(prev_ramp.values,dsize=(width,length), interpolation=cv2.INTER_LINEAR) #or INTER_CUBIC ?
        # if gacos is to be applied, need to remove its effect first here:
        if 'gacos' in ifg_ml.variables:
            resized = resized - ifg_ml['gacos'].values
        ifg_ml['toremove'] = ifg_ml['toremove'] + resized
        # apply the correction - wrapping it and conjugate in complex realm
        mag = np.ones(resized.shape)
        correction = wrap2phase(resized)
        cpx_corr = magpha2RI_array(mag, correction)
        #a trick to apply the correction only to non-nan values.. 'other'
        #da = da.where(xrda.isnull(), other=da.values * np.conjugate(cpx_corr))
        ifg_ml['cpx'].values = ifg_ml['cpx'].values * np.conjugate(cpx_corr)
        resized = ''
        if pre_detrend:
            print('no need to detrend if prev_ramp is here, cancelling to avoid extra noise')
            pre_detrend = False
    ifg_ml['pha'] = ifg_ml.cpx
    ifg_ml['pha'].values = np.angle(ifg_ml.cpx)
    #if type(ml_weights) == type(None):
    #    #have the orig coh here:
    #    ifg_ml['coh'] = ifg_ml['pha']
    #    ifg_ml.coh.values = np.abs(ifg_ml.cpx)
    #have gacos removed first, prior to doing height corr:
    if 'gacos' in ifg_ml.variables:
        pha_no_gacos = wrap2phase(ifg_ml['pha'] - ifg_ml['gacos'])
        #if np.nanstd(pha_no_gacos) >= np.nanstd(ifg_ml.pha.values):
        if get_fft_std(pha_no_gacos) >= get_fft_std(ifg_ml['pha'].values):
            print('GACOS correction would increase overall phase std - dropping')
            #ifg_ml = ifg_ml.drop('gacos')
        else:
            ifg_ml['pha'].values = pha_no_gacos #wrap2phase(ifg_ml['pha'] - ifg_ml['gacos'])
            ifg_ml['pha'] = ifg_ml['pha'].where(ifg_ml.mask>0)
            ifg_ml['toremove'] = ifg_ml['toremove'] + ifg_ml['gacos']
            ifg_ml['gacos'] = ifg_ml['gacos'].where(ifg_ml.mask>0)
            #ok, return coh, phase back to cpx
            cpxa = magpha2RI_array(ifg_ml.coh.values, ifg_ml.pha.values)
            ifg_ml['cpx'].values = cpxa
    #
    if pre_detrend:
        ifg_ml['cpx'], correction = detrend_ifg_xr(ifg_ml['cpx'], isphase=False, return_correction = True)
        ifg_ml['pha'].values = np.angle(ifg_ml.cpx)
        #ifg_ml['pha'] = ifg_ml['pha'].where(ifg_ml.mask>0)
        ifg_ml['toremove'] = ifg_ml['toremove'] + correction
    
    # just mask it
    ifg_ml['pha'] = ifg_ml['pha'].where(ifg_ml.mask>0)
    ifg_ml['coh'] = ifg_ml['coh'].where(ifg_ml.mask>0)
    
    if 'hgt' in ifg.variables:
        ifg_ml['hgt'] = ifg_ml['hgt'].where(ifg_ml.mask>0)
    # perform Gaussian filtering
    #ifg_ml = filter_ifg_ml(ifg_ml)
    #now fix the correlation with heights:
    if hgtcorr:
        #ifg_ml['toremove'] = ifg_ml['toremove'] + 
        # dounw=False may be faster!
        print('calculating correlation with DEM')
        try:
            toremove_hgt = remove_height_corr(ifg_ml, tmpdir = tmpdir, dounw = True, nonlinear=False)
            pha_no_hgt = wrap2phase(ifg_ml['pha'] - toremove_hgt)
            if np.nanstd(pha_no_hgt) >= np.nanstd(ifg_ml.pha.values):
                print('but the correction would increase overall phase std - dropping')
                hgtcorr = False
            else:
                ifg_ml['toremove'] = ifg_ml['toremove'] + toremove_hgt
                # need to remove hgt only here, as the 'toremove' was already removed before..
                ifg_ml['pha'].values = pha_no_hgt
                ifg_ml['pha'] = ifg_ml['pha'].where(ifg_ml.mask>0)
                #ok, return coh, phase back to cpx
                cpxa = magpha2RI_array(ifg_ml.coh.values, ifg_ml.pha.values)
                ifg_ml['cpx'].values = cpxa
                if pre_detrend:
                    ifg_ml['cpx'], correction = detrend_ifg_xr(ifg_ml['cpx'], isphase=False, return_correction = True)
                    ifg_ml['pha'].values = np.angle(ifg_ml.cpx)
                    ifg_ml['pha'] = ifg_ml['pha'].where(ifg_ml.mask>0)
                    ifg_ml['toremove'] = ifg_ml['toremove'] + correction
        except:
            print('some error trying correlate with DEM. continuing without it')
            # ifg_ml['pha'].plot()
    #else:
    #    ifg_ml['toremove'] = 0*ifg_ml['pha']
    #here we remove height correlation - gacos was already removed from pha
    #ifg_ml['pha'].values = wrap2phase(ifg_ml['pha'] - ifg_ml['toremove'])
    #ifg_ml['pha'] = ifg_ml['pha'].where(ifg_ml.mask>0)
    # maybe not the best, but have gacos correction inside the toremove variable
    ifg_ml['toremove'] = ifg_ml['toremove'].where(ifg_ml.mask>0)
    # oh, ok, also masking pixels with small number of pre-multilooked points (landmasked)
    if 'pxcount' in ifg_ml.data_vars:
        # try setting to something like 90 if ML10
        if not thres_pxcount:
            #thres_pxcount = int(round((ml**2)/2))
            thres_pxcount = int(round((ml**2)*4/5))
        ifg_ml['mask'] = ifg_ml.mask * (ifg_ml.pxcount >= thres_pxcount)
        ifg_ml['pha'] = ifg_ml['pha'].where(ifg_ml.mask>0)
        ifg_ml['cpx'] = ifg_ml['cpx'].where(ifg_ml.mask>0)
        ifg_ml['coh'] = ifg_ml['coh'].where(ifg_ml.mask>0)
        #ifg_ml['gacos'] = ifg_ml['gacos'].where(ifg_ml.mask>0)
    #print('finally, filter using (adapted) gauss filter')
    if ml > 2:
        calc_coh_from_delta = True
    else:
        # that part takes ages and it is not that big improvement..
        calc_coh_from_delta = False
    ifg_ml = filter_ifg_ml(ifg_ml, calc_coh_from_delta = calc_coh_from_delta)
    if not keep_coh_debug:
        ifg_ml['coh'] = ifg_ml['orig_coh']
        ifg_ml['coh'] = ifg_ml['coh'].where(ifg_ml.mask>0)
    return ifg_ml


def load_tif2xr(tif):
    #xrpha = xr.open_rasterio(tif)
    xrpha = rioxarray.open_rasterio(tif)
    xrpha = xrpha.squeeze('band')
    xrpha = xrpha.drop('band')
    xrpha = xrpha.rename({'x': 'lon','y': 'lat'})
    return xrpha


def deramp_unw(xrda, dim=['lat','lon']):
    da = xrda.fillna(0).copy(deep=True)
    dt = xr.apply_ufunc(
                    _detrend_2d_ufunc,
                    da,
                    input_core_dims=[dim],
                    output_core_dims=[dim],
                    output_dtypes=[da.dtype],
                    vectorize=True,
                    dask="parallelized",
                )
    return dt


def _detrend_2d_ufunc(arr):
    assert arr.ndim == 2
    N = arr.shape
    col0 = np.ones(N[0] * N[1])
    col1 = np.repeat(np.arange(N[0]), N[1]) + 1
    col2 = np.tile(np.arange(N[1]), N[0]) + 1
    G = np.stack([col0, col1, col2]).transpose()
    d_obs = np.reshape(arr, (N[0] * N[1], 1))
    m_est = np.dot(np.dot(spl.inv(np.dot(G.T, G)), G.T), d_obs)
    d_est = np.dot(G, m_est)
    linear_fit = np.reshape(d_est, N)
    return arr - linear_fit


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
    
    


def detrend_ifg_xr(xrda, isphase=True, return_correction = False, maxfringes = 4):
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
    if (abs(numfringesx) > maxfringes) or (abs(numfringesy) > maxfringes):
        print('too many fringes identified, probably wrong. cancelling detrend')
        if return_correction:
            return da, 0
        else:
            return da
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



    
def filter_ifg_ml(ifg_ml, calc_coh_from_delta = False):  #, rotate = False):
    # this will use only PHA and would filter it
    #normalise mag
    tempar_mag1 = np.ones_like(ifg_ml.pha)
    ifg_ml['cpx'].values = magpha2RI_array(tempar_mag1, ifg_ml.pha.values)
    print('filter using (adapted) gauss filter')
    ifg_ml['gauss_cpx'] = filter_cpx_gauss(ifg_ml, sigma = 1, trunc = 2)
    '''
    # tested - rotate doesn't do any difference - the gauss filter is ok in this way!
    if rotate:
        rotpha = wrap2phase(ifg_ml.pha.values + np.pi/4)
        ifg_ml['cpx'].values = magpha2RI_array(tempar_mag1, rotpha)
        ifg_ml['gauss_cpx2'] = filter_cpx_gauss(ifg_ml, sigma = 1, trunc = 2)
        ifg_ml['gauss_pha2'] = 0*ifg_ml['pha']
        ifg_ml['gauss_pha2'].values = np.angle(ifg_ml.gauss_cpx2.values)
        ifg_ml['gauss_pha2'].values = wrap2phase(ifg_ml.gauss_pha2.values - np.pi/4)
        #return orig cpx back
        ifg_ml['cpx'].values = magpha2RI_array(tempar_mag1, ifg_ml.pha.values)
        # and now somehow get the max of both filters - think about how - abs/pha or just.. sin/cos maxes?
        r=np.real(ifg_ml['gauss_cpx2'].values)
        i=np.imag(ifg_ml['gauss_cpx2'].values)
        # ..............
    '''
    ifg_ml['gauss_pha'] = 0*ifg_ml['pha']
    ifg_ml['gauss_pha'].values = np.angle(ifg_ml.gauss_cpx.values)
    #use magnitude after filtering as coherence
    ifg_ml['gauss_coh'] = 0*ifg_ml['pha']
    # that is great but has problems at maxima of cos or sin
    ifg_ml['gauss_coh'].values = np.abs(ifg_ml.gauss_cpx.values)
    if calc_coh_from_delta:
        # if using the coh from the phase residuals (based on variance), it 
        delta = np.angle(np.exp(1j*(ifg_ml['gauss_pha'] - ifg_ml['pha'])))
        #aaa = coh_from_phadiff(delta, winsize=5)
        print('calculating coh from phase diff')
        phacoh = coh_from_phadiff(delta, winsize=3)
        # this should be much better:
        #ifg_ml['gauss_coh'].values = 1-np.abs(delta)/np.pi
        #but it is not actually. so i am using max of gauss_coh and phacoh:
        ifg_ml['gauss_coh'].values=np.maximum(ifg_ml['gauss_coh'].values, phacoh)
    return ifg_ml


def block_unwrap(da_pha, da_coh):
    return da_unw


def block_correlate(da1, da2):
    return dacorr


def wrap2phase(A):
    return np.angle(np.exp(1j*A))


import time

# 2022-04-04 - ok, this should now work after the updated, properly!
def cascade_unwrap(frame, pair, downtoml = 1, procdir = os.getcwd(), only10 = True, smooth = False, hgtcorr = True, outtif = None, subtract_gacos = False, cliparea_geo = None, dolocal = False):
    '''
    only10 = only 1 previous ramp, scaled 10x to the downtoml, 
    20220506 - i turned 'smooth' off as it caused large errors!!! not sure why, need to revisit the filtering!!!!
    '''
    print('performing cascade unwrapping')
    starttime = time.time()
    if only10:
        #ifg_ml10 = process_ifg(frame, pair, procdir = procdir, ml = 10, fillby = 'gauss', defomax = 0.5, add_resid = False, smooth = True)
        #ifg_ml1 = process_ifg(frame, pair, procdir = procdir, ml = 1, fillby = 'gauss', prevest = ifg_ml10['unw'], defomax = 0.6, add_resid = True, smooth=True, outtif=outtif)
        # orig:
        #ifg_ml10 = process_ifg(frame, pair, procdir = procdir, ml = 10, fillby = 'gauss', defomax = 0.5, add_resid = False, hgtcorr = hgtcorr, rampit=True)
        # 01/2022: updating parameters:
        ifg_ml10 = process_ifg(frame, pair, procdir = procdir, ml = 10*downtoml, fillby = 'gauss', defomax = 0.3, thres = 0.4, add_resid = False, hgtcorr = hgtcorr, rampit=True, dolocal = dolocal)
        if downtoml == 1:
            # avoiding gauss proc, as seems heavy for memory
            ifg_ml = process_ifg(frame, pair, procdir = procdir, ml = downtoml, fillby = 'nearest', smooth = False, prev_ramp = ifg_ml10['unw'], defomax = 0.3, thres = 0.3, add_resid = True, hgtcorr = False, outtif=outtif, subtract_gacos = subtract_gacos, cliparea_geo = cliparea_geo,  dolocal = dolocal)
        else:
            ifg_ml = process_ifg(frame, pair, procdir = procdir, ml = downtoml, fillby = 'gauss', prev_ramp = ifg_ml10['unw'], defomax = 0.3, thres = 0.3, add_resid = True, hgtcorr = False, outtif=outtif, subtract_gacos = subtract_gacos, cliparea_geo = cliparea_geo,  dolocal = dolocal, smooth=smooth)
    else:
        ifg_mlc = process_ifg(frame, pair, procdir = procdir, ml = 20, fillby = 'gauss', defomax = 0.5, add_resid = False, hgtcorr = hgtcorr, rampit=True,  dolocal = dolocal)
        for i in [10, 5, 3]:
            if downtoml < i:
                ifg_mla = process_ifg(frame, pair, procdir = procdir, ml = i, fillby = 'gauss', prev_ramp = ifg_mlc['unw'], defomax = 0.5, add_resid = False, hgtcorr = hgtcorr, rampit=True,  dolocal = dolocal)
                ifg_mlc = ifg_mla.copy(deep=True)
        if downtoml == 1:
            # avoiding gauss proc, as seems heavy for memory
            ifg_ml = process_ifg(frame, pair, procdir = procdir, ml = downtoml, fillby = 'nearest', smooth = False, prev_ramp = ifg_mlc['unw'], defomax = 0.3, thres = 0.3, add_resid = True, hgtcorr = False, outtif=outtif, subtract_gacos = subtract_gacos,  dolocal = dolocal)
        else:
            ifg_ml = process_ifg(frame, pair, procdir = procdir, ml = downtoml, fillby = 'gauss', prev_ramp = ifg_mlc['unw'], thres = 0.3, defomax = 0.4, add_resid = True, hgtcorr = False, outtif=outtif, subtract_gacos = subtract_gacos, cliparea_geo = cliparea_geo,  dolocal = dolocal)
        #ifg_ml10 = process_ifg(frame, pair, procdir = procdir, ml = 10, fillby = 'gauss', prevest = ifg_ml20['unw'], defomax = 0.5, add_resid = False, rampit=True)
        ##elapsed_time10 = time.time()-starttime
        #ifg_ml5 = process_ifg(frame, pair, procdir = procdir, ml = 5, fillby = 'gauss', prevest = ifg_ml10['unw'], defomax = 0.6, add_resid = False, rampit=True)
        #ifg_ml3 = process_ifg(frame, pair, procdir = procdir, ml = 3, fillby = 'gauss', prevest = ifg_ml5['unw'], defomax = 0.6, add_resid = False, rampit=True)
        #ifg_ml1 = process_ifg(frame, pair, procdir = procdir, ml = 1, fillby = 'gauss', prevest = ifg_ml3['unw'], defomax = 0.6, add_resid = True, outtif=outtif)
    elapsed_time = time.time()-starttime
    hour = int(elapsed_time/3600)
    minite = int(np.mod((elapsed_time/60),60))
    sec = int(np.mod(elapsed_time,60))
    print("\nTotal elapsed time: {0:02}h {1:02}m {2:02}s".format(hour,minite,sec))
    return ifg_ml

#print('preparing nc and png')
#bin2nc(outunwbin, mask_full, outunwnc, dtype = np.float32)
#create_preview(outunwnc, ftype = 'unwrapped')




def make_avg_amp(mlitiflist, hgtxr):
    avgamp = hgtxr*0
    nopixels = avgamp.copy()
    for ampf in mlitiflist:
        try:
            amp = io.read_geotiff(ampf) #/ 255
        except:
            print('error reading '+ampf)
            continue
        amp[np.isnan(amp)] = 0
        nopixels.values[amp>0] += 1
        avgamp = avgamp + amp
    avgamp = avgamp/nopixels
    return avgamp


def make_std_amp(mlitiflist, avgamp):
    ddof = 1
    avgvar = avgamp*0
    nopixels = avgvar.copy()
    for ampf in mlitiflist:
        try:
            amp = io.read_geotiff(ampf)
        except:
            print('error reading '+ampf)
            continue
        amp[np.isnan(amp)] = 0
        nopixels.values[amp>0] += 1
        avgvar = avgvar + (amp - avgamp)**2
    # correct for ddof
    sumpx = nopixels - ddof
    #sumpx[sumpx<1] = np.nan
    sumpx.values[sumpx<1] = np.nan
    avgstd = np.sqrt(avgvar/sumpx)
    return avgstd

'''
# for unlimited time dimension netcdf

import netCDF4

def _expand_variable(nc_variable, data, expanding_dim, nc_shape, added_size):
    # For time deltas, we must ensure that we use the same encoding as
    # what was previously stored.
    # We likely need to do this as well for variables that had custom
    # econdings too
    if hasattr(nc_variable, 'calendar'):
        data.encoding = {
            'units': nc_variable.units,
            'calendar': nc_variable.calendar,
        }
    data_encoded = xr.conventions.encode_cf_variable(data) # , name=name)
    left_slices = data.dims.index(expanding_dim)
    right_slices = data.ndim - left_slices - 1
    nc_slice   = (slice(None),) * left_slices + (slice(nc_shape, nc_shape + added_size),) + (slice(None),) * (right_slices)
    nc_variable[nc_slice] = data_encoded.data



def append_to_netcdf(filename, ds_to_append, unlimited_dims):
    if isinstance(unlimited_dims, str):
        unlimited_dims = [unlimited_dims]
    if len(unlimited_dims) != 1:
        # TODO: change this so it can support multiple expanding dims
        raise ValueError(
            "We only support one unlimited dim for now, "
            f"got {len(unlimited_dims)}.")
    unlimited_dims = list(set(unlimited_dims))
    expanding_dim = unlimited_dims[0]
    with netCDF4.Dataset(filename, mode='a') as nc:
        nc_dims = set(nc.dimensions.keys())
        nc_coord = nc[expanding_dim]
        nc_shape = len(nc_coord)
        added_size = len(ds_to_append[expanding_dim])
        variables, attrs = xr.conventions.encode_dataset_coordinates(ds_to_append)
        for name, data in variables.items():
            if expanding_dim not in data.dims:
                # Nothing to do, data assumed to the identical
                continue
            nc_variable = nc[name]
            _expand_variable(nc_variable, data, expanding_dim, nc_shape, added_size)


# TODO - make median coh, but better - RAM-effective - using dask
def make_median_coh_memory(group):
    #medcoh = hgtxr*0
    firstpass = True
    #thirddim = 'pair'
    for pair,row in group.iterrows():
        cohf = row['cohf']
        try:
            data = xr.open_rasterio(cohf)
        except:
            print('error reading coh of pair '+pair)
            continue
        if firstpass:
            medcoh_xr = data.copy(deep=True).squeeze('band').drop('band')
            coh_np = data.values
            firstpass = False
        else:
            try:
                coh_np = np.append(coh_np, data.values, axis = 0)
            except:
                print('error in pair '+pair+' - maybe shape mismatch?')
    medcoh_xr.values = np.nanmedian(coh_np, axis=0)
    number_of_cohs = coh_np.shape[0]
    return medcoh_xr, number_of_cohs
'''


'''
import dask
def make_median_coh(group, num_workers = 1, outtif = None):
    #outnc is only a temp file
    outnc = 'tmpcoh.nc'
    varname = 'coherence'
    if os.path.exists(outnc):
        os.remove(outnc)
    #medcoh = hgtxr*0
    firstpass = True
    thirddim = 'pair'
    for pair,row in group.iterrows():
        cohf = row['cohf']
        try:
            data = xr.open_rasterio(cohf)
            data = data.rename({'x': 'lon','y': 'lat','band': thirddim})
            data[thirddim] = [pair]
            #sort lons and lats - some tifs may have it inverted?
            data = data.sortby([thirddim,'lon','lat'])
        except:
            print('error reading coh of pair '+pair)
            continue
        if firstpass:
            #coh_xr = ds.copy(deep=True) #.squeeze('band').drop('band')
            #coh_np = data.values
            ds = xr.Dataset({varname:data})
            coordsys = data.crs.split('=')[1]
            ds = ds.rio.set_spatial_dims(x_dim="lon", y_dim="lat", inplace=True)
            ds = ds.rio.write_crs(coordsys, inplace=True)
            ds.to_netcdf(outnc, mode='w', unlimited_dims=[thirddim])
            firstpass = False
        else:
            try:
                #fast fix having sometimes wrong lon/lats (still same dimensions)
                data['lon'] = ds.lon
                data['lat'] = ds.lat
                dsA = xr.Dataset({varname:data})
                append_to_netcdf(outnc, dsA, unlimited_dims=thirddim)
                #coh_np = np.append(coh_np, data.values, axis = 0)
            except:
                print('error in pair '+pair+' - maybe shape mismatch?')
    if not os.path.exists(outnc):
        print('some error, cancelling for frame '+frame)
        return False
    print('temporary chunked netcdf ready, now dask-calculate median')
    dsA = ''  #clean from RAM
    data = ''
    # this will load it to dask - accessible through 'data' later
    medcoh_xr = xr.open_dataset(outnc, chunks={"lon":10,"lat":10})
    medcoh_dask = dask.array.nanmedian(medcoh_xr.coherence.data, axis=0)
    #xr.apply_ufunc(
    #    np.nanmedian,
    #    medcoh_xr.coherence,
    #    input_core_dims=[[dim]], # [dim]],
    #    dask="parallelized",
    #    kwargs={"axis": 0})
    #    output_dtypes=[np.uint8],
    #)
    medcoh = medcoh_dask.compute(num_workers=num_workers)  # scheduler='processes' # output is numpy ndarray
    #save to netcdf:
    ds = ds.squeeze(thirddim).drop(thirddim)
    da = ds['coherence']
    da.values = medcoh.astype(np.uint8)
    da.attrs['NUMBER_OF_INPUT_FILES'] = len(medcoh_xr[thirddim])
    if outtif:
        da = da.rio.set_spatial_dims(x_dim="lon", y_dim="lat", inplace=True)
        da = da.rio.write_crs(coordsys, inplace=True)
        da.rio.to_raster(outtif, compress='deflate')
    return da

'''
def make_avg_coh(group, hgtxr):
    avgcoh = hgtxr*0
    nopixels = avgcoh.copy()
    for cohf in group['cohf']:
        try:
            coh = io.read_geotiff(cohf) / 255
        except:
            print('error reading '+cohf)
            continue
        coh[np.isnan(coh)] = 0
        nopixels.values[coh>0] += 1
        avgcoh = avgcoh + coh
    avgcoh = avgcoh/nopixels
    return avgcoh


def make_std_coh(group, avgcoh):
    ddof = 1
    avgvar = avgcoh*0
    nopixels = avgvar.copy()
    for cohf in group['cohf']:
        try:
            coh = io.read_geotiff(cohf) / 255
        except:
            print('error reading '+cohf)
            continue
        coh[np.isnan(coh)] = 0
        nopixels.values[coh>0] += 1
        avgvar = avgvar + (coh - avgcoh)**2
    # correct for ddof
    sumpx = nopixels - ddof
    #sumpx[sumpx<1] = np.nan
    sumpx.values[sumpx<1] = np.nan
    avgstd = np.sqrt(avgvar/sumpx)
    return avgstd


def get_date_matrix(pairs):
    date_matrix = pd.DataFrame(pairs, columns=['pair'])
    date_matrix['date1'] = pd.to_datetime(date_matrix.pair.str[:8], format='%Y%m%d')
    date_matrix['date2'] = pd.to_datetime(date_matrix.pair.str[9:], format='%Y%m%d')
    date_matrix['btemp'] = date_matrix.date2 - date_matrix.date1
    date_matrix = date_matrix.set_index('pair')
    return date_matrix


def build_amp_avg_std(frame, return_ampstab = False):
    track=str(int(frame[:3]))
    epochsdir = os.path.join(os.environ['LiCSAR_public'], track, frame, 'epochs')
    hgtfile = os.path.join(os.environ['LiCSAR_public'], track, frame, 'metadata', frame+'.geo.hgt.tif')
    hgtxr = xr.open_rasterio(os.path.join(hgtfile))
    hgtxr = hgtxr.squeeze('band').drop('band')
    mlitiflist = fc.get_epochs(frame, return_mli_tifs = True)
    print('generating amp average')
    ampavg = make_avg_amp(mlitiflist, hgtxr)
    print('generating amp std')
    ampstd = make_std_amp(mlitiflist, ampavg)
    if return_ampstab:
        ampstab = 1 - ampstd/ampavg
        ampstab.values[ampstab<=0] = 0.00001
        return ampstab
    else:
        return ampavg, ampstd


def build_coh_avg_std(frame, ifgdir = None, days = 'all', monthly = False, outnopx = False, do_std = True, do_tif = False):
    track=str(int(frame[:3]))
    if not ifgdir:
        ifgdir = os.path.join(os.environ['LiCSAR_public'], track, frame, 'interferograms')
    try:
        pairs = get_ifgdates(ifgdir)
    except:
        print('error, dropping frame '+frame)
        pairs = None
        return False
    pairsall = get_date_matrix(pairs)
    pairsall['pair'] = list(pairsall.index)
    pairsall['cohf'] = ifgdir + '/' + pairsall.pair + '/' + pairsall.pair + '.geo.cc.tif'
    if days != 'all':
        pairs = pairsall[pairsall.btemp == str(days)+' days']
    #pairs['datetocheck'] = pairs.date1 + pd.Timedelta('6 days')
    #pairs['month'] = pairs['datetocheck'].dt.month
    hgtfile = os.path.join(os.environ['LiCSAR_public'], track, frame, 'metadata', frame+'.geo.hgt.tif')
    hgtxr = xr.open_rasterio(os.path.join(hgtfile))
    hgtxr = hgtxr.squeeze('band').drop('band')
    if not monthly:
        print('generating coh average')
        if days == 'all':
            print('fast written solution, will output xr Dataset of avg cohs')
            coherences = xr.Dataset()
            #ndays = pairsall.btemp.dt.days.unique()
            #ndays.sort()
            daygroups = [(6,6), (12,12), (18,24), (30,42), (48,72), (78,96), (102,156), (162,200), (201,300), (301,400)]
            for daygroup in daygroups:
                print('preparing coh avgs for btemp between {0} and {1} days'.format(daygroup[0],daygroup[1]))
                pairs = pairsall[np.isin(pairsall.btemp.dt.days, np.arange(daygroup[0], daygroup[1]+1))]
                #now, this will select in periods between 1st March and 1st September
                for setX in [([3,4,5,6,7,8], 'summer'), ([9,10,11,12,1,2], 'winter')]:
                    setX_time = setX[0]
                    setX_label = setX[1]
                    center_date = pairs.date1+(pairs.date2-pairs.date1)/2
                    setA = pairs[np.isin(center_date.dt.month,setX_time)]
                    if len(setA) > 1:
                        cohavgA = make_avg_coh(setA, hgtxr)
                        cohavgA.attrs['number of input cohs'] = len(setA)
                        nameA = 'mean coh {0}-{1} days {2}'.format(daygroup[0], daygroup[1], setX_label)
                        coherences[nameA] = cohavgA
                # cohavgB.plot(vmin=0,vmax=0.9)
                # plt.show()
                #for ddays in days:
                # 
            return coherences
        else:
            cohavg = make_avg_coh(pairs, hgtxr)
            if do_std:
                print('generating coh std')
                cohstd = make_std_coh(pairs, cohavg)
            else:
                cohstd = 0
    else:
        pairs['month'] = pairs['date1'].dt.month
        for i, group in pairs.groupby('month'):
            print('preparing month '+str(i)+' from {} coh maps'.format(str(len(group))))
            out = make_avg_coh(group, hgtxr)
            outtif = frame+'.avg_coh.month'+str(i)+'.tif'
            print('exporting to '+outtif)
            out.rio.write_crs("epsg:4326", inplace=True).rio.to_raster(outtif)
    if do_tif:
        outtif = os.path.join(os.environ['LiCSAR_public'], track, frame, 'metadata', frame+'.geo.meancoh.'+str(days)+'.tif')
        nopx = len(pairs)
        cohavg = (cohavg*255).astype(np.uint8)
        cohavg.attrs['NUMBER_OF_INPUT_FILES'] = nopx
        cohavg = cohavg.rename({'x': 'lon','y': 'lat'}).sortby(['lon','lat'])
        export_xr2tif(cohavg, outtif, debug = False)
    if outnopx:
        nopx = len(pairs)
        return cohavg, cohstd, nopx
    else:
        return cohavg, cohstd



'''
# once in /work/scratch-pw/earmla/LiCSAR_temp/batchdir/142D_09148_131313
coherences = build_coh_avg_std(frame, ifgdir = 'GEOC', days = 'all')
coherences.to_netcdf('cohnewmeans.nc')
# to plot them
i=0
plt.close()
for x in coherences:
    i=i+1
    label = str(x)+': from {} files'.format(coherences[x].attrs['number of input cohs'])
    name = 'meancoh_{}.png'.format(str(i))
    coherences[x].plot(vmin=0,vmax=0.9)
    plt.title(label)
    print(label)
    plt.savefig(name, format='png')
    plt.close()
    #plt.show()
    
'''


'''
def build_median_coh(frame, days = 12, num_workers = 8):
    track=str(int(frame[:3]))
    ifgdir = os.path.join(os.environ['LiCSAR_public'], track, frame, 'interferograms')
    try:
        pairs = get_ifgdates(ifgdir)
    except:
        print('error, dropping frame '+frame)
        pairs = None
        return False
    pairs = get_date_matrix(pairs)
    pairs = pairs[pairs.btemp == str(days)+' days']
    pairs['pair'] = list(pairs.index)
    #pairs['datetocheck'] = pairs.date1 + pd.Timedelta('6 days')
    #pairs['month'] = pairs['datetocheck'].dt.month
    pairs['cohf'] = ifgdir + '/' + pairs.pair + '/' + pairs.pair + '.geo.cc.tif'
    outtif = os.path.join(os.environ['LiCSAR_public'], track, frame, 'metadata', frame+'.geo.medcoh.'+str(days)+'.tif')
    if os.path.exists(outtif):
        print('already exists: '+outtif)
    print('generating median of {} days coh for frame '.format(str(days))+frame)
    da = make_median_coh(pairs, num_workers = num_workers)
    export_xr2tif(da, outtif, debug = False)



def build_median_coh_frames(frames, dayset = [6, 12, 18], avg = True):
    lenf = len(frames)
    i = 0
    for frame in frames['frameID']:
        i = i+1
        print('[{0}/{1}] processing frame {2}'.format(str(i),str(lenf),frame))
        for days in dayset:
            try:
                if avg:
                    dmm = build_coh_avg_std(frame, days = days, monthly = False, outnopx = False, do_std = False, do_tif = True)
                    dmm = ''
                else:
                    build_median_coh(frame, days = days)
            except:
                print('error processing '+frame)

'''

'''

tr=134
orb='D'
days = 12

import LiCSAR_lib.unwrp_multiscale as unw
#for tr in range(134,140):
for tr in range(100,176):
    #for tr in range(150,160):
    #for tr in range(160,176):
    frames=unw.fc.lq.sqlout2list(unw.fc.lq.get_frames_in_orbit(tr,orb))
    for frame in frames:
        try:
            dmm = unw.build_coh_avg_std(frame, days = days, monthly = False, outnopx = False, do_std = False, do_tif = True)
        except:
            print('error in '+frame)


    unw.build_median_coh_frames(frames, dayset, avg = True)


def build_median_coh_all_frames(dayset = [6, 12, 18], avg = True):
    ascf, descf = unw.fc.get_all_frames()
    6,
    12
    , 
    [24] 
    dayset = [36,48]
    frames = descf
    for frames in [ascf, descf]:
        unw.build_median_coh_frames(frames, dayset, avg = True)
        
'''



'''
how to make multiprocessing? this way:
from pathos.multiprocessing import ProcessingPool as Pool
import multiprocessing as multi

p = Pool(n_para_gap)
_result = np.array(p.map(count_gaps_wrapper, range(n_para_gap)),
                                   dtype=object)

p.close()



q = multi.get_context('fork')

p = q.Pool(n_para_gap)
                _result = np.array(p.map(count_gaps_wrapper, range(n_para_gap)),
                                   dtype=object)
                p.close()


'''
