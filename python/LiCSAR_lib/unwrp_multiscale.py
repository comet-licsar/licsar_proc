################################################################################
#Imports
################################################################################
import numpy as np
import os
import subprocess
import xarray as xr
import cv2
from scipy import interpolate
from scipy import ndimage
import time
import matplotlib.pyplot as plt
xr.set_options(keep_attrs=True)

from scipy.ndimage import gaussian_filter
from scipy.ndimage import generic_filter
from scipy import stats

from astropy.convolution import Gaussian2DKernel, interpolate_replace_nans, convolve_fft

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



def main_unwrap(cpxbin, cohbin, maskbin, outunwbin, width, est = None):
    print('processing by snaphu')
    prefix = os.path.dirname(cpxbin)+'/'
    snaphuconffile = make_snaphu_conf(prefix)
    extracmd = ''
    if est:
        extracmd = "-e {}".format(est)
    starttime = time.time()
    snaphucmd = 'snaphu -f {0} -M {1} -o {2} -c {3} {4} {5} {6}'.format(snaphuconffile, maskbin, outunwbin, cohbin, extracmd, cpxbin, str(width))
    runcmd(snaphucmd,True)
    elapsed_time = time.time()-starttime
    hour = int(elapsed_time/3600)
    minite = int(np.mod((elapsed_time/60),60))
    sec = int(np.mod(elapsed_time,60))
    print("\nElapsed time: {0:02}h {1:02}m {2:02}s".format(hour,minite,sec))


def create_preview(ncfile, ftype = 'unwrapped'):
    if 'wrapped' not in ftype:
        print('wrong ftype')
        return False
    tosource = os.path.join(os.environ['LiCSARpath'],'lib','LiCSAR_bash_lib.sh')
    command = 'create_preview_'+ftype
    os.system('source {0}; {1} {2}'.format(tosource, command, ncfile))





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


def make_snaphu_conf(sdir):
    snaphuconf = ('STATCOSTMODE  DEFO\n',
            'INFILEFORMAT  COMPLEX_DATA\n',
            'CORRFILEFORMAT  FLOAT_DATA\n'
            'OUTFILEFORMAT FLOAT_DATA\n',
            'ESTFILEFORMAT FLOAT_DATA\n',
            'RMTMPTILE TRUE\n')
    snaphuconffile = os.path.join(sdir,'snaphu.conf')
    with open(snaphuconffile,'w') as f:
        for l in snaphuconf:
            f.write(l)
    return snaphuconffile


def process_ifg(frame, pair, procdir = os.getcwd(), ml = 10, fillby = 'gauss'):
    pubdir = os.environ['LiCSAR_public']
    geoframedir = os.path.join(pubdir,str(int(frame[:3])),frame)
    geoifgdir = os.path.join(geoframedir,'interferograms',pair)
    tmpdir = os.path.join(procdir,pair,'temp_'+str(ml))
    if not os.path.exists(os.path.join(procdir,pair)):
        os.mkdir(os.path.join(procdir,pair))
    if not os.path.exists(tmpdir):
        os.mkdir(tmpdir)
    #orig files
    ifg_pha_file = os.path.join(geoifgdir,pair+'.geo.diff_pha.tif')
    coh_file = os.path.join(geoifgdir,pair+'.geo.cc.tif')
    landmask_file = os.path.join(geoframedir,'metadata',frame+'.geo.landmask.tif')
    #orig file ncs
    origphanc = os.path.join(tmpdir,'origpha.nc')
    origcohnc = os.path.join(tmpdir,'origcoh.nc')
    landmasknc = os.path.join(procdir,'landmask.nc')
    # prepare input dataset and load
    convert_tif2nc(ifg_pha_file,origphanc)
    #this one is tricky - would cause issues if will try to multiprocess..
    if not os.path.exists(landmasknc):
        convert_tif2nc(landmask_file,landmasknc)
    convert_tif2nc(coh_file,origcohnc, coh = True)
    #use orig ifg and coh
    inpha = xr.open_dataset(origphanc)
    incoh = xr.open_dataset(origcohnc)
    inmask = incoh.copy(deep=True)
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
    ifg['coh'].values = incoh.z.values
    ifg['mask'].values = inmask.z.values
    inpha = ''
    incoh = ''
    inmask = ''
    # now doing multilooking, using coh as mag...
    #make complex from coh and pha
    ifg['cpx'] = ifg.coh.copy()
    ifg['cpx'].values = magpha2RI_array(ifg.coh.values, ifg.pha.values)
    ifg_ml = multilook_by_gauss(ifg, ml)
    #mask it ... this worked ok here:
    #thres = 0.5
    #gmmm.... it works better if i use the 'trick on gauss of normalised ifg....'
    thres = 0.35
    ifg_ml['mask_coh'] = ifg_ml['mask'].where(ifg_ml.gauss_coh > thres).fillna(0)
    #try without modal filter..
    #ifg_ml = filter_mask_modal(ifg_ml, 'mask_coh', 'mask_coh', 8)
    print('interpolate coh-based masked areas of gauss pha')
    tofillpha = ifg_ml.pha.fillna(0).where(ifg_ml.mask_coh.where(ifg_ml.mask == 1).fillna(1) == 1)
    if fillby == 'gauss':
        #trying astropy approach now:
        print('filling through Gaussian kernel')
        kernel = Gaussian2DKernel(x_stddev=2)
        # create a "fixed" image with NaNs replaced by interpolated values
        ifg_ml['pha'].values = interpolate_replace_nans(tofillpha.values, kernel)
    else:
        print('filling by nearest neighbours')
        ifg_ml['pha'].values = interpolate_nans(tofillpha.values, method='nearest')
        #ifg_ml['gauss_pha'] = ifg_ml['gauss_pha'].fillna(0)
    #exporting for snaphu
    #normalise mag from the final pha
    tempar_mag1 = np.ones_like(ifg_ml.pha)
    cpxarr = magpha2RI_array(tempar_mag1, ifg_ml.pha.fillna(0).values)
    ifg_ml['gauss_cpx'].values = cpxarr
    print('unwrapping by snaphu')
    binmask= os.path.join(tmpdir,'gaussmask.bin')
    bincoh = os.path.join(tmpdir,'gausscoh.bin')
    binR = os.path.join(tmpdir,'gaussR.bin')
    binI = os.path.join(tmpdir,'gaussI.bin')
    binCPX = os.path.join(tmpdir,'cpxgaussifg.bin')
    outunwbin = os.path.join(tmpdir,'gaussunwrapped.bin')
    #print('exporting to bin files')
    ifg_ml.mask_coh.fillna(0).values.astype(np.byte).tofile(binmask)
    ifg_ml.gauss_coh.fillna(0).values.astype(np.float32).tofile(bincoh)
    #ifg_ml.coh.fillna(0).values.astype(np.float32).tofile(bincoh)
    r = np.real(ifg_ml.gauss_cpx).astype(np.float32).fillna(0).values.tofile(binR)
    i = np.imag(ifg_ml.gauss_cpx).astype(np.float32).fillna(0).values.tofile(binI)
    RI2cpx(binR, binI, binCPX)
    width = len(ifg_ml.lon)
    #unwrapping itself
    main_unwrap(binCPX, bincoh, binmask, outunwbin, width)
    print('importing snaphu result to ifg_ml')
    binfile = outunwbin
    #toxr = ifg_ml
    daname = 'unw'
    dtype = np.float32
    unw1 = np.fromfile(binfile,dtype=dtype)
    unw1 = unw1.reshape(ifg_ml.pha.shape)
    #unw1 = np.flip(unw1,axis=0)
    ifg_ml[daname] = ifg_ml['pha'].copy()
    ifg_ml[daname].values = unw1
    #ok, so the gauss-based coh mask is not the best to do... so exporting 'all pixels'
    ifg_ml[daname] = ifg_ml[daname]*ifg_ml['mask']
    #ifg_ml[daname] = ifg_ml[daname]*ifg_ml['mask_coh']
    ifg_ml[daname].values = ifg_ml[daname].values - np.nanmedian(ifg_ml[daname].where(ifg_ml.mask_coh>0).values)
    ifg_ml[daname] = ifg_ml[daname]*ifg_ml.mask_coh
    return ifg_ml


def process_frame(frame, ml = 10):
    #the best to run in directory named by the frame id
    pubdir = os.environ['LiCSAR_public']
    geoframedir = os.path.join(pubdir,str(int(frame[:3])),frame)
    geoifgdir = os.path.join(geoframedir,'interferograms')
    for pair in os.listdir(geoifgdir):
        if os.path.exists(os.path.join(geoifgdir, pair, pair+'.geo.diff_pha.tif')):
            if not os.path.exists(os.path.join(pair,pair+'.unw')):
                print('processing pair '+pair)
                ifg_ml = process_ifg(frame, pair, procdir = os.getcwd(), ml = ml, fillby = 'gauss')
                #ifg_ml.unw.values.tofile(pair+'/'+pair+'.unw')
                #np.flipud(ifg_ml.unw.fillna(0).values).tofile(pair+'/'+pair+'.unw')
                np.flipud(ifg_ml.unw.where(ifg_ml.mask_coh > 0).values).tofile(pair+'/'+pair+'.unw')
                #or use gauss here?
                np.flipud((ifg_ml.coh.where(ifg_ml.mask > 0)*255).astype(np.byte).fillna(0).values).tofile(pair+'/'+pair+'.cc')
                width = len(ifg_ml.lon)
                create_preview_bin(pair+'/'+pair+'.unw', width, ftype = 'unw')
                os.system('rm '+pair+'/'+pair+'.unw.ras')
                os.system('rm -r '+pair+'/'+'temp_'+str(ml))


def multilook_by_gauss(ifg, ml = 10):
    kernel = Gaussian2DKernel(x_stddev=ml)
    resultR = convolve_fft(np.real(ifg['cpx'].values), kernel)
    resultI = convolve_fft(np.imag(ifg['cpx'].values), kernel)
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


def multilook_normalised(ifg, ml = 10):
    #landmask it and multilook it
    bagcpx = ifg[['cpx']].where(ifg.mask>0).coarsen({'lat': ml, 'lon': ml}, boundary='trim')
    ifg_ml = bagcpx.sum() / bagcpx.count()
    ifg_ml['pha'] = ifg_ml.cpx.copy()
    ifg_ml['pha'].values = np.angle(ifg_ml.cpx)
    #have the orig coh here:
    ifg_ml['coh'] = ifg_ml['pha']
    ifg_ml.coh.values = np.abs(ifg_ml.cpx)
    #downsample mask
    ifg_ml['mask'] = ifg.mask.coarsen({'lat': ml, 'lon': ml}, boundary='trim').max()
    #normalise mag
    tempar_mag1 = np.ones_like(ifg_ml.pha)
    cpxarr = magpha2RI_array(tempar_mag1, ifg_ml.pha.values)
    ifg_ml['cpx'].values = cpxarr
    print('filter using (adapted) gauss filter')
    ifg_ml['gauss_cpx'] = filter_cpx_gauss(ifg_ml, sigma = 1, trunc = 2)
    ifg_ml['gauss_pha'] = ifg_ml.pha
    ifg_ml['gauss_pha'].values = np.angle(ifg_ml.gauss_cpx.values)
    #use magnitude after filtering as coherence
    ifg_ml['gauss_coh'] = ifg_ml.pha
    ifg_ml['gauss_coh'].values = np.abs(ifg_ml.gauss_cpx.values)
    return ifg_ml

