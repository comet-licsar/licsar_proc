#!/usr/bin/env python3

################################################################################
# LiCSAR Unwrapper
# by Milan Lazecky, 2021-2022, University of Leeds
#
# version: 1.0.0 (2022-06-03)
#
# A tool to unwrap LiCSAR (or any other) interferogram, starting from geotiffs
# Mandatory inputs: geotiffs of phase, coherence
# Optional inputs: geotiffs with GACOS corrections, DEM, landmask (automatically found for LiCSAR)
#
# Pre-requisities: snaphu
# Optional requisites: GMT, cpxfiddle (doris), ImageMagick
#
################################################################################
#Imports
################################################################################
import os, glob
import shutil
import subprocess

import xarray as xr
xr.set_options(keep_attrs=True)
import rioxarray
from osgeo import gdal

import numpy as np
import pandas as pd

from scipy import interpolate
from scipy import ndimage
from scipy.ndimage import gaussian_filter
from scipy.ndimage import generic_filter
from scipy import stats
import scipy.signal as sps
import scipy.linalg as spl
from astropy.convolution import Gaussian2DKernel, interpolate_replace_nans, convolve_fft
from sklearn.linear_model import HuberRegressor

import time
import matplotlib.pyplot as plt


# avoid cv2 in ipynb
def in_ipynb():
    try:
        cfg = get_ipython().config 
        return True
    except NameError:
        return False


if not in_ipynb():
    # some extra imports, used by additional functions
    try:
        import cv2
    except:
        print('cv2 not loaded - cascade will not work')
else:
    print('at JASMIN notebook service, cv2 does not load - cascade will not work')


try:
    import dask.array as da
except:
    print('dask not loaded - hgt correlation will not work')


try:
    import framecare as fc
except:
    print('framecare not loaded - amplitude stability will not work')


try:
    from LiCSAR_lib.LiCSAR_misc import *
except:
    print('licsar misc not loaded')


try:
    import LiCSBAS_io_lib as io
    from LiCSBAS_tools_lib import *
except:
    print('licsbas not loaded - the amplitude/coherence average/stability will fail')


################################################################################
# Main functions to perform the unwrapping
################################################################################

def cascade_unwrap(frame, pair, downtoml = 1, procdir = os.getcwd(),
                   only10 = True, smooth = False, thres=0.3, hgtcorr = True, 
                   outtif = None, cliparea_geo = None, subtract_gacos = False, dolocal = False):
    """Main function to unwrap a geocoded LiCSAR interferogram using a cascade approach.
    
    Args:
        frame (string): LiCSAR frame ID
        pair (string): identifier of interferometric pair, e.g. '20200120_20200201'
        downtoml (int): target multilook factor (default: 1, no extra multilooking)
        procdir (string): path to processing directory
        only10 (boolean): switch to use only 1 previous ramp, scaled 10x to the downtoml, instead of few cascades
        smooth (boolean): switch to use extra Gaussian filtering for 2-pass unwrapping
        thres (float): threshold between 0-1 for gaussian-based coherence-like measure (spatial phase consistence?); higher number - more is masked prior to unwrapping
        
        hgtcorr (boolean): switch to perform correction for height-phase correlation.
        outtif (string or None): path to geotiff file to export result to.
        cliparea_geo (string or None): use GMT/LiCSBAS string to identify area to clip, in geo-coordinates, as ``'lon1/lon2/lat1/lat2'``
        subtract_gacos (boolean): switch whether to return the interferograms with GACOS being subtracted (by default, GACOS is used only to support unwrapping and would be added back)
        dolocal (boolean): switch to use local directory to find interferograms, rather than search for LiCSAR_public directory in JASMIN
    
    Returns:
        xarray.Dataset: unwrapped multilooked interferogram with additional layers
    """
    print('performing cascade unwrapping')
    starttime = time.time()
    if only10:
        # 01/2022: updating parameters:
        ifg_ml10 = process_ifg(frame, pair, procdir = procdir, ml = 10*downtoml, fillby = 'gauss', defomax = 0.3, thres = 0.4, add_resid = False, hgtcorr = hgtcorr, rampit=True, dolocal = dolocal, smooth=True)
        if downtoml == 1:
            # avoiding gauss proc, as seems heavy for memory
            ifg_ml = process_ifg(frame, pair, procdir = procdir, ml = downtoml, fillby = 'nearest', smooth = False, prev_ramp = ifg_ml10['unw'], defomax = 0.3, thres = thres, add_resid = True, hgtcorr = False, outtif=outtif, subtract_gacos = subtract_gacos, cliparea_geo = cliparea_geo,  dolocal = dolocal)
        else:
            ifg_ml = process_ifg(frame, pair, procdir = procdir, ml = downtoml, fillby = 'gauss', prev_ramp = ifg_ml10['unw'], defomax = 0.3, thres = thres, add_resid = True, hgtcorr = False, outtif=outtif, subtract_gacos = subtract_gacos, cliparea_geo = cliparea_geo,  dolocal = dolocal, smooth=smooth)
    else:
        ifg_mlc = process_ifg(frame, pair, procdir = procdir, ml = 20, fillby = 'gauss', defomax = 0.5, add_resid = False, hgtcorr = hgtcorr, rampit=True,  dolocal = dolocal)
        for i in [10, 5, 3]:
            if downtoml < i:
                ifg_mla = process_ifg(frame, pair, procdir = procdir, ml = i, fillby = 'gauss', prev_ramp = ifg_mlc['unw'], defomax = 0.5, add_resid = False, hgtcorr = hgtcorr, rampit=True,  dolocal = dolocal)
                ifg_mlc = ifg_mla.copy(deep=True)
        if downtoml == 1:
            # avoiding gauss proc, as seems heavy for memory
            ifg_ml = process_ifg(frame, pair, procdir = procdir, ml = downtoml, fillby = 'nearest', smooth = False, prev_ramp = ifg_mlc['unw'], defomax = 0.3, thres = thres, add_resid = True, hgtcorr = False, outtif=outtif, subtract_gacos = subtract_gacos,  dolocal = dolocal)
        else:
            ifg_ml = process_ifg(frame, pair, procdir = procdir, ml = downtoml, fillby = 'gauss', prev_ramp = ifg_mlc['unw'], thres = thres, defomax = 0.4, add_resid = True, hgtcorr = False, outtif=outtif, subtract_gacos = subtract_gacos, cliparea_geo = cliparea_geo,  dolocal = dolocal)
    elapsed_time = time.time()-starttime
    hour = int(elapsed_time/3600)
    minite = int(np.mod((elapsed_time/60),60))
    sec = int(np.mod(elapsed_time,60))
    print("\nTotal elapsed time: {0:02}h {1:02}m {2:02}s".format(hour,minite,sec))
    return ifg_ml


def load_tif(frame,pair,dtype='unw',cliparea_geo=None):
    pubdir = os.environ['LiCSAR_public']
    geoframedir = os.path.join(pubdir,str(int(frame[:3])),frame)
    geoifgdir = os.path.join(geoframedir,'interferograms',pair)
    infile = os.path.join(geoifgdir,pair+'.geo.'+dtype+'.tif')
    return  load_tif2xr(infile,cliparea_geo=cliparea_geo)


def get_resolution(ifg, in_m=True):
    """Gets resolution of the xr.dataset (or dataarray), either in metres or degrees
    """
    resdeg = (np.abs(ifg.lat[1]-ifg.lat[0])+np.abs(ifg.lon[1]-ifg.lon[0]))/2
    if in_m:
        latres = 111.32 * np.cos(np.radians(ifg.lat.mean())) * 1000 # in m
        return float(latres * resdeg)
    else:
        return float(resdeg)


def load_ifg(frame, pair, unw=True, dolocal=False):
    pubdir = os.environ['LiCSAR_public']
    geoframedir = os.path.join(pubdir,str(int(frame[:3])),frame)
    if dolocal:
        geoifgdir = os.path.join('GEOC',pair)
    else:
        geoifgdir = os.path.join(geoframedir,'interferograms',pair)
    #orig files
    # will use only the filtered ifgs now..
    ifg_pha_file = os.path.join(geoifgdir,pair+'.geo.diff_pha.tif')
    coh_file = os.path.join(geoifgdir,pair+'.geo.cc.tif')
    landmask_file = os.path.join(geoframedir,'metadata',frame+'.geo.landmask.tif')
    hgtfile = os.path.join(geoframedir,'metadata',frame+'.geo.hgt.tif')
    # load the files
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
    # create datacube
    ifg = xr.Dataset()
    ifg['pha'] = inpha
    ifg['coh'] = ifg['pha']
    ifg['coh'].values = incoh.values
    ifg['mask'] = ifg['pha']
    ifg['mask'].values = inmask.values
    # just to clean from memory
    inpha=''
    incoh=''
    # to load orig unw_file
    if unw:
        unw_file = os.path.join(geoifgdir,pair+'.geo.unw.tif')
        incoh = load_tif2xr(unw_file)
        ifg['unw'] = ifg['pha']
        ifg['unw'].values = inunw.values
    ifg['mask_extent'] = ifg['pha'].where(ifg['pha'] == 0).fillna(1)
    # including hgt anyway - would be useful later
    if os.path.exists(hgtfile):
        try:
            inhgt = load_tif2xr(hgtfile)
            ifg['hgt'] = ifg['pha']
            ifg['hgt'].values = inhgt.values
        except:
            print('ERROR in importing heights!')
            hgtcorr = False
    return ifg


def process_ifg(frame, pair, procdir = os.getcwd(), 
        ml = 10, fillby = 'gauss', thres = 0.35, smooth = False, defomax = 0.3,
        hgtcorr = False, gacoscorr = True, pre_detrend = True,
        cliparea_geo = None, outtif = None, prevest = None, prev_ramp = None,
        coh2var = True, add_resid = True,  rampit=False, subtract_gacos = False, dolocal = False,
        cohratio = None, keep_coh_debug = True):
    """Main function to unwrap a geocoded LiCSAR interferogram. Works on JASMIN (but can be easily adapted for local use)

    Args:
        frame (string): LiCSAR frame ID
        pair (string): identifier of interferometric pair, e.g. ``'20200120_20200201'``
        procdir (string): path to processing directory
        ml (int): multilooking factor used to reduce the interferogram in lon/lat
        fillby (string): algorithm to fill gaps. use one of values: ``'gauss'``, ``'nearest'``, ``'none'`` (where ``'none'`` would only fill NaNs by zeroes)
        thres (float): threshold between 0-1 for gaussian-based coherence-like measure (spatial phase consistence?); higher number - more is masked prior to unwrapping
        smooth (boolean): switch to use extra Gaussian filtering for 2-pass unwrapping
        defomax (float): parameter to snaphu for maximum deformation in rad per 2pi cycle (DEFOMAX_CYCLE)
        
        hgtcorr (boolean): switch to perform correction for height-phase correlation
        gacoscorr (boolean): switch to apply GACOS corrections (if detected)
        pre_detrend (boolean): switch to apply detrending on wrapped phase to support unwrapping
        
        cliparea_geo (string or None): use GMT/LiCSBAS string to identify area to clip, in geo-coordinates, as ``'lon1/lon2/lat1/lat2'``
        outtif (string or None): path to geotiff file to export result to.
        prevest (xarray.DataArray or None): a previous rough estimate to be used by snaphu as the ESTFILE
        prev_ramp (xarray.DataArray or None): a previous estimate or a ramp that will be removed prior to unwrapping (and added back)
        
        coh2var (boolean): convert coherence to variance for weighting. could be useful, but need to change from squared, something to try...
        add_resid (boolean): switch to add back residuals from spatially filtered unwrapping (makes sense if smooth is ON)
        rampit (boolean): perform an extra strong gaussian filter to get a very rough unwrapping result. basically a longwave signal ramp. used by cascade approach
        subtract_gacos (boolean): switch whether to return the interferograms with GACOS being subtracted (by default, GACOS is used only to support unwrapping and would be added back)
        dolocal (boolean): switch to use local directory to find interferograms, rather than search for LiCSAR_public directory in JASMIN
        
        cohratio (xr.DataArray): coherence ratio (or another array) to be used for weighting the phase instead of the original coherence
        keep_coh_debug (boolean): only in combination with use_coh_stab - whether or not to keep original (downsampled) ifg coherence after using the coh_stab to weight the phase during multilooking
    
    Returns:
        xarray.Dataset: unwrapped multilooked interferogram with additional layers
    """
    try:
        ifg = load_ifg(frame, pair, dolocal=dolocal)
    except:
        print('error in loading data')
        return False
    
    # prepare tmp dir structure
    tmpdir = os.path.join(procdir,pair,'temp_'+str(ml))
    tmpgendir = os.path.join(procdir,pair,'temp_gen')
    if not os.path.exists(os.path.join(procdir,pair)):
        os.mkdir(os.path.join(procdir,pair))
    if not os.path.exists(tmpdir):
        os.mkdir(tmpdir)
    if not os.path.exists(tmpgendir):
        os.mkdir(tmpgendir)
    
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
        ifg['gacos'] = ifg.pha
        ifg['gacos'].values = ingacos.values
    else:
        gacoscorr = False
    # masking by coherence if we do not use multilooking - here the coherence corresponds to reality
    if ml == 1:
        cohthres = 0.15
        ifg['mask'] = ifg['mask'] * ifg['mask'].where(ifg['coh'] < cohthres).fillna(1)
    if hgtcorr:
        if not 'hgt' in ifg:
            print('ERROR in importing heights!')
            hgtcorr = False
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
    # now crop if needed:
    if cliparea_geo:
        minclipx, maxclipx, minclipy, maxclipy = cliparea_geo.split('/')
        minclipx, maxclipx, minclipy, maxclipy = float(minclipx), float(maxclipx), float(minclipy), float(maxclipy)
        if minclipy > maxclipy:
            print('you switched min max in crop coordinates (latitude). fixing')
            tmpcl = minclipy
            minclipy=maxclipy
            maxclipy=tmpcl
        if minclipx > maxclipx:
            print('you switched min max in crop coordinates (longitude). fixing')
            tmpcl = minclipx
            minclipx=maxclipx
            maxclipx=tmpcl
        # now will clip it - lat is opposite-sorted, so need to slice from max to min in y ... plus 10 pixels on all sides
        resdeg = get_resolution(ifg, in_m=False)
        ifg = ifg.sel(lon=slice(minclipx-10*resdeg, maxclipx+10*resdeg), lat=slice(maxclipy+10*resdeg, minclipy-10*resdeg))
        # not the best here, as pixels might get slightly shifted, but perhaps not that big deal (anyway prev_ramp is 'blurred')
        if not type(prev_ramp) == type(None):
            prev_ramp = prev_ramp.sel(lon=slice(minclipx-10*resdeg, maxclipx+10*resdeg), lat=slice(maxclipy+10*resdeg, minclipy-10*resdeg))
    #WARNING - ONLY THIS FUNCTION HAS GACOS INCLUDED NOW! (and heights fix!!!)
    ifg_ml = multilook_normalised(ifg, ml, tmpdir = tmpdir, hgtcorr = hgtcorr, pre_detrend = pre_detrend, prev_ramp = prev_ramp, keep_coh_debug = keep_coh_debug)
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
    lenthres = 7000 # m
    # resolution of orig ifg is expected 0.1 km
    mlres = get_resolution(ifg_ml, in_m=True)
    pixels = int(round(lenthres/mlres))
    #origres = 0.1
    #resdeg = np.abs(ifg.lat[1]-ifg.lat[0])
    #latres = 111.32 * np.cos(np.radians(ifg_ml.lat.mean()))
    #origres = float(latres * resdeg) # in km
    #
    #pixels = int(round(lenthres/origres/ml))
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
        #ifg_ml = filter_ifg_ml(ifg_ml)
        # 2022/07: adding strong filter, say radius 1 km
        ifg_ml = filter_ifg_ml(ifg_ml, radius = 1000)
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


def process_frame(frame, ml = 10, thres = 0.35, smooth = False, cascade=False, 
            hgtcorr = True, gacoscorr = True,
            cliparea_geo = None, pairsetfile = None, 
            export_to_tif = False, subtract_gacos = False,
            nproc = 1, dolocal = False,
            use_amp_stab = False, use_coh_stab = False, keep_coh_debug = True):
    """Main function to process whole LiCSAR frame (i.e. unwrap all available interferograms within the frame). Works only at JASMIN.

    Args:
        frame (string): LiCSAR frame ID
        ml (int): multilooking factor used to reduce the interferogram in lon/lat
        thres (float): threshold between 0-1 for gaussian-based coherence-like measure (spatial phase consistence?); higher number - more is masked prior to unwrapping
        smooth (boolean): switch to use extra Gaussian filtering for 2-pass unwrapping
        cascade (boolean): switch to perform cascade unwrapping
        
        hgtcorr (boolean): switch to perform correction for height-phase correlation
        gacoscorr (boolean): switch to apply GACOS corrections (if detected)
        
        cliparea_geo (string or None): use GMT/LiCSBAS string to identify area to clip, in geo-coordinates: e.g. ``'lon1/lon2/lat1/lat2'``
        pairsetfile (string or None): path to file containing list of pairs to unwrap
        export_to_tif (boolean): switch to export unwrapped data to geotiffs (default: False, generate only binaries, as used by LiCSBAS)
        subtract_gacos (boolean): switch whether to return the interferograms with GACOS being subtracted (by default, GACOS is used only to support unwrapping and would be added back)
        nproc (int): use multiprocessing (one core per interferogram), not well tested, uses pathos
        dolocal (boolean): switch to use local directory to find interferograms, rather than search for LiCSAR_public directory in JASMIN
        
        use_amp_stab (boolean): apply amplitude stability index instead of coherence-per-interferogram for unwrapping
        use_coh_stab (boolean): apply (experimental) coherence stability index. not recommended (seems not logical to me) - worth investigating though (maybe helps against loop closure errors)
        keep_coh_debug (boolean): only in combination with use_coh_stab - whether or not to keep original (downsampled) ifg coherence after using the coh_stab to weight the phase during multilooking
    
    Returns:
        xarray.Dataset: multilooked interferogram with additional layers
    """
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
    # functions for multiprocessing
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
                        defomax = 0.3
                        ifg_ml = process_ifg(frame, pair, procdir = os.getcwd(), ml = ml, hgtcorr = hgtcorr, fillby = 'gauss', 
                                 thres = thres, defomax = defomax, add_resid = True, outtif = outtif, cohratio = cohratio, smooth = smooth,
                                 keep_coh_debug = keep_coh_debug, gacoscorr = gacoscorr, cliparea_geo = cliparea_geo,
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
                            ifg_ml = cascade_unwrap(frame, pair, downtoml = ml, procdir = os.getcwd(), outtif = outtif, subtract_gacos = subtract_gacos, smooth = smooth, cliparea_geo = cliparea_geo, dolocal = dolocal)
                        else:
                            defomax = 0.3
                            ifg_ml = process_ifg(frame, pair, procdir = os.getcwd(), ml = ml, hgtcorr = hgtcorr, fillby = 'gauss', 
                                     thres = thres, defomax = defomax, add_resid = True, outtif = outtif, cohratio = cohratio, 
                                     keep_coh_debug = keep_coh_debug, gacoscorr = gacoscorr, cliparea_geo = cliparea_geo,
                                     subtract_gacos = subtract_gacos, dolocal = dolocal, smooth = smooth)
                        (ifg_ml.unw.where(ifg_ml.mask_full > 0).values).astype(np.float32).tofile(pair+'/'+pair+'.unw')
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
                        os.system('rm -r '+pair+'/'+'temp_*')
                    except:
                        print('ERROR processing of pair '+pair)
                        os.system('rm -r '+pair)
                if not os.path.exists(os.path.join(pair,pair+'.unw')):
                    print('some error occured and the unw was not processed')
                    os.system('rm -r '+pair)



def multilook_normalised(ifg, ml = 10, tmpdir = os.getcwd(), hgtcorr = True, pre_detrend = True, prev_ramp = None, thres_pxcount = None, keep_coh_debug = True):
    """Multilooking function that does much more.
    
    This function is normally called by process_ifg. It would use coherence as weights to multilook interferometric phase, and downsample other layers if available to a final datacube.
    It will apply mask, including based on number of valid pixels in the multilooking window.
    It will also apply Gaussian filter, mainly to get the Gaussian-based coherence-like measure (used no matter if smooth is ON).

    Args:
        ifg (xarray.Dataset): xarray dataset containing interferogram layers, mainly cpx for complex numbers interferogram
        ml (int): multilooking factor used to reduce the interferogram in both x/y or lon/lat
        tmpdir (string): path to temporary directory
        hgtcorr (boolean): switch to perform correction for height-phase correlation
        pre_detrend (boolean): switch to perform detrending of phase
        prev_ramp (xarray.DataArray): a previous (ramp) estimate. it can be of different dimensions as it would get interpolated
        thres_pxcount (int or None): by default, we nullify multilooked pixel that has less than 4/5 non-nan input values. You may change this, e.g. if ml=10, apply thres_pxcount=90 for keeping only pixel with over 9/10 values
        keep_coh_debug (boolean): for experiments, this would keep the original interferogram coherence instead of use average coherence or amplitude stability etc. for weighting

    Returns:
        xarray.Dataset: multilooked interferogram with additional layers
    """
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
        ''' removing this check, because we want to FORCE-apply GACOS.. otherwise we get loop closure errors...
        pha_no_gacos = wrap2phase(ifg_ml['pha'] - ifg_ml['gacos'])
        #if np.nanstd(pha_no_gacos) >= np.nanstd(ifg_ml.pha.values):
        if get_fft_std(pha_no_gacos) >= get_fft_std(ifg_ml['pha'].values):
            print('GACOS correction would increase overall phase std - dropping')
            #ifg_ml = ifg_ml.drop('gacos')
        else:
            ifg_ml['pha'].values = pha_no_gacos
        '''
        ifg_ml['pha'].values = wrap2phase(ifg_ml['pha'] - ifg_ml['gacos'])
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
    #
    # just mask it
    ifg_ml['pha'] = ifg_ml['pha'].where(ifg_ml.mask>0)
    ifg_ml['coh'] = ifg_ml['coh'].where(ifg_ml.mask>0)
    #
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
            cpx_no_hgt = magpha2RI_array(np.abs(ifg_ml.cpx.values), pha_no_hgt)
            #if np.nanstd(pha_no_hgt) >= np.nanstd(ifg_ml.pha.values):
            #    print('but the correction would increase overall phase std - dropping')
            #    hgtcorr = False
            if np.nanstd(cpx_no_hgt) >= np.nanstd(ifg_ml.cpx.values):
                print('but the correction would increase overall complex std - dropping')
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
    # maybe not the best, but have gacos correction inside the toremove variable
    ifg_ml['toremove'] = ifg_ml['toremove'].where(ifg_ml.mask>0)
    # oh, ok, also masking pixels with small number of pre-multilooked points (landmasked)
    if 'pxcount' in ifg_ml.data_vars:
        # try setting to something like 90 if ML10
        if not thres_pxcount:
            #thres_pxcount = int(round((ml**2)/2))
            # if not set, we will auto-set it to mask multilooked pixels with less than 80% input (non-nan) pixels
            thres_pxcount = int(round((ml**2)*4/5))
        ifg_ml['mask'] = ifg_ml.mask * (ifg_ml.pxcount >= thres_pxcount)
        ifg_ml['pha'] = ifg_ml['pha'].where(ifg_ml.mask>0)
        ifg_ml['cpx'] = ifg_ml['cpx'].where(ifg_ml.mask>0)
        ifg_ml['coh'] = ifg_ml['coh'].where(ifg_ml.mask>0)
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


################################################################################
# Helping functions
################################################################################

def interpolate_nans(array, method='cubic'):
    """Interpolation of NaN values in a grid

    Args:
        array (np.array): numpy array with nans to interpolate
        method (string): interpolation method for griddata function, e.g. cubic

    Returns:
        np.array: interpolated grid
    """
    array = np.ma.masked_invalid(array)
    x = np.arange(0, array.shape[1])
    y = np.arange(0, array.shape[0])
    xx, yy = np.meshgrid(x, y)
    x1 = xx[~array.mask]
    y1 = yy[~array.mask]
    newarr = array[~array.mask]
    GD1 = interpolate.griddata((x1, y1), newarr.ravel(),(xx, yy),method=method)
    GD1 = np.array(GD1)
    return GD1


def runcmd(cmd, printcmd = True):
    """Runs command through os.system

    Args:
        cmd (string): command to run
        printcmd (boolean): if True, will do verbose
    """
    if printcmd:
        print(cmd)
        rc = os.system(cmd)
    else:
        #with nostdout():
        rc = os.system(cmd+' >/dev/null 2>/dev/null')
    if not rc == 0:
        print('WARNING - command did not exit as OK')


def magpha2RI_array(mag, pha):
    """Converts arrays of magnitude and phase to complex number array (real and imaginary)

    Args:
        mag (np.array): numpy array with magnitude values
        pha (np.array): numpy array with phase values

    Returns:
        np.array: complex number array
    """
    R = np.cos(pha) * mag
    I = np.sin(pha) * mag
    out = R + 1j*I
    return out


def coh_from_phadiff(phadiff, winsize = 3):
    """Calculates coherence based on variance of interferogram, computed in window with given size

    Args:
        phadiff (np.array): interferogram
        winsize (int): window size

    Returns:
        np.array: coherence based on the variance
    """
    variance = ndimage.generic_filter(phadiff, np.var, size=winsize)
    outcoh = 1/np.sqrt(1+winsize*winsize*variance)
    return outcoh


def filter_cpx_gauss(ifg_ml, sigma = 2, trunc = 4):
    """Gaussian-based spatial filter on complex numbers (interferogram)

    Args:
        ifg_ml (xr.Dataset): xarray dataset of the interferogram, must contain 'cpx' dataarray
        sigma (int): sigma parameter to gaussian filter
        trunc (int): trunc parameter to gaussian filter

    Returns:
        xr.Dataarray: filtered complex numbers dataarray
    """
    # tried with R, I separately --- EXACT same result as if using cpx numbers...
    #R = np.real(ifg_ml.cpx.values)
    #I = np.imag(ifg_ml.cpx.values)
    #
    #gR = filter_nan_gaussian_conserving(R, sigma=sigma, trunc=trunc)
    #gI = filter_nan_gaussian_conserving(I, sigma=sigma, trunc=trunc)
    #
    #gauss_cpx = gR + 1j*gI
    gauss_cpx = filter_nan_gaussian_conserving(ifg_ml.cpx.values, sigma=sigma, trunc=trunc)
    gauss_xr = ifg_ml['cpx'].copy(deep=True)
    #
    gauss_xr.values = gauss_cpx
    return gauss_xr


def filter_nan_gaussian_conserving(arr, sigma, trunc):
    """Apply a gaussian filter to an array with nans.
    
    based on:
    https://stackoverflow.com/questions/18697532/gaussian-filtering-a-image-with-nan-in-python
    
    Intensity is only shifted between not-nan pixels and is hence conserved.
    The intensity redistribution with respect to each single point
    is done by the weights of available pixels according
    to a gaussian distribution.
    All nans in arr, stay nans in gauss.
    
    Args:
        arr (np.array): array of real numbers to filter by Gaussian kernel
        sigma (int): sigma parameter to gaussian filter
        trunc (int): trunc parameter to gaussian filter

    Returns:
        np.array: filtered real numbers array
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


def create_preview_bin(binfile, width, ftype = 'unw'):
    """Use of cpxfiddle to create simple preview PNG rasters from binary files

    Args:
        binfile (string): path to the binary file to generate preview from
        width (int): width of the binary file
        ftype (string): filetype. Supported types: unw, pha, coh, mag

    Returns:
        string: filename of the generated preview png file
    """
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


try:
    def resize_bin(inbin, inwid, inlen, outbin, outwid, outlen, dtype = np.byte, intertype = cv2.INTER_NEAREST):
        """Use of cv2 to interpolate/resize binary file to new dimensions

        Args:
            inbin (string): path to the binary file
            inwid (int): width of the input binary file
            inlen (int): length of the input binary file
            outbin (string): path to the output binary file
            outwid (int): target width of the output binary file
            outlen (int): target length of the output binary file
            dtype (string): data type of the binary, e.g. np.byte
            intertype (string): type of interpolation for cv2, e.g. cv2.INTER_NEAREST, cv2.INTER_LINEAR, cv2.INTER_CUBIC
        """
        a = np.fromfile(inbin, dtype=dtype).reshape(inlen, inwid)
        #use cv2.INTER_CUBIC for upsample np.float32 data ...
        #a = a - np.nanmedian(a)
        out = cv2.resize(a,dsize=(outwid,outlen), interpolation=intertype)
        out.astype(dtype).tofile(outbin)
        return
except:
    print('error loading resize_bin function - cascade will not work (install cv2)')



def RI2cpx(R, I, cpxfile):
    """Convert real and imaginary binary files to a complex number binary file. Obsolete function.
    
    Args:
        R (string): path to the binary file with real values
        I (string): path to the binary file with imaginary values
        cpxfile (string): path to the binary file for complex output
    """
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


def remove_islands(npa, pixelsno = 50):
    '''Removes isolated clusters of pixels from numpy array npa having less than pixelsno pixels.
    
    Args:
        npa (np.array): (unwrapped) interferogram with NaNs
        pixelsno (int): minimum number of pixels in isolated clusters (connected components)
    
    Returns:
        np.array: array after removing islands
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


def main_unwrap(cpxbin, cohbin, maskbin = None, outunwbin = 'unwrapped.bin', width = 0, est = None, bin_pre_remove = None, defomax = 0.6, printout = True):
    '''Main function to perform unwrapping with snaphu.
    
    Args:
        cpxbin (string): path to cpxfloat32 binary interferogram to unwrap
        cohbin (string): path to float32 binary for coherence
        maskbin (string or None): path to mask binary
        outunwbin (string): path to output unwrapped binary
        width (int): width of binary raster
        est (string or None): path to coarse estimate binary (float32)
        bin_pre_remove (string or None): path to float32 binary to remove from est, prior to unwrapping
        defomax (float): max defo cycles
        printout (boolean): controls verbosity of text output
    '''
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
    """Creates preview of interferogram (wrapped or unwrapped) - works only with licsar_proc
    
    Args:
        infile (string): path to input tif to generate preview
        ftype (string): type of the input file. can be: 'wrapped', 'unwrapped'
    """
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


def make_snaphu_conf(sdir, defomax = 1.2):
    """Creates snaphu configuration file
    
    Args:
        sdir (string): directory where to generate snaphu.conf
        defomax (float): DEFOMAX parameter to snaphu
    
    Returns:
        string: path to generated snaphu.conf
    """
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
    """Creates GACOS correction for the interferogram. works only at JASMIN
    
    Args:
        frame (string): frame ID
        pair (string): pair ID (e.g. '20201001_20201201')
    
    Returns:
        string: path to generated GACOS correction, or False
    """
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
    '''Removes height-correlated signal
    
     first, correlate ```ifg_ml['pha']``` and ```ifg_ml['hgt']``` in blocks
      - better to keep dounw=True that unwraps each block by snaphu. but it can be slow
     get coefficient for correction of correlating areas
     interpolate the coefficient throughout whole raster
     multiply by hgt = 'to_remove'
     
     Args:
        ifg_ml (xarray.Dataset): input dataset
        corr_thres (float): threshold of correlation within window to keep
        tmpdir (string): path to temp directory
        dounw (boolean): whether to unwrap the small windows, or try correlate with wrapped phase
        nonlinear (boolean): pass nonlinear parameter to correct_hgt
    
    Returns:
        xarray.Dataarray: dataarray of height correlation corrections
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
    # first unwrap the block if conditions are ok
    toret = None
    try:
        hgt=hgt.values
    except:
        print('ok, hgt was already np')
    if hgt[hgt != 0].size == 0:
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
            unwr = unwrap_np(cpx.values, coh.values, defomax = 0.3, tmpdir=tmpdir)
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


def unwrap_xr(ifg, mask=True, defomax = 0.3, tmpdir=os.getcwd()):
    """Quite direct unwrapping of the xarray Dataset of ifg
    
    Args:
        ifg (xarray.Dataset): ifg dataset
        mask (boolean): whether to use mask
        defomax (float): DEFOMAX to snaphu
        tmpdir (string): temp dir
    
    Returns:
        xarray.Dataset: ifg dataset now with unwrapped result
    """
    coh = ifg.coh.values
    cpx = ifg.cpx.fillna(0).astype(np.complex64).values
    if mask:
        mask = ifg.mask.fillna(0).values
    unw = unwrap_np(cpx, coh, mask = mask, defomax = defomax, tmpdir=tmpdir, deltemp=False)
    ifg['unw']=ifg.pha.copy(deep=True)
    ifg['unw'].values = unw
    return ifg


def unwrap_np(cpx, coh, defomax = 0.3, tmpdir=os.getcwd(), mask = False, deltemp = True):
    '''unwraps given numpy array
    
    Args:
        cpx (numpy.ndarray): array of complex interferogram
        coh (numpy.ndarray): array of coherence
        defomax (float): DEFOMAX to snaphu
        tmpdir (string): temp dir
        mask (boolean): whether to try use binary mask (if exists)
        deltemp (boolean): clean temp dir after processing
    
    Returns:
        numpy.ndarray: unwrapped array
    '''
    if not os.path.exists(tmpdir):
        os.mkdir(tmpdir)
    if mask:
        try:
            binmask= os.path.join(tmpdir,'mask.bin')
            mask.astype(np.byte).tofile(binmask)
        except:
            binmask=None
    else:
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


def export_xr2tif(xrda, tif, lonlat = True, debug = True, dogdal = True):
    """Exports xarray dataarray to a geotiff
    
     Args:
        xrda (xarray.Dataarray): dataarray to export
        tif (string): path to output tif file
        lonlat (boolean): are the dimensions named as lon, lat?
        debug (boolean): just load it as float32
        dogdal (boolean): after exporting, perform gdalwarp (fix for potential issues in output geotiff)
    """
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


def load_tif2xr(tif, cliparea_geo=None, tolonlat=True):
    """loads geotiff to xarray.DataArray
    
    Args:
        tif (string): path to geotiff
        cliparea_geo (string): use GMT/LiCSBAS string to identify area to clip, in geo-coordinates, as ``'lon1/lon2/lat1/lat2'``
        tolonlat (boolean): if True, return as lon lat coordinates
    
    Returns:
        xr.DataArray: loaded contents
    """
    xrpha = rioxarray.open_rasterio(tif)
    xrpha = xrpha.squeeze('band')
    xrpha = xrpha.drop('band')
    
    if cliparea_geo:
        minclipx, maxclipx, minclipy, maxclipy = cliparea_geo.split('/')
        minclipx, maxclipx, minclipy, maxclipy = float(minclipx), float(maxclipx), float(minclipy), float(maxclipy)
        if minclipy > maxclipy:
            print('you switched min max in crop coordinates (latitude). fixing')
            tmpcl = minclipy
            minclipy=maxclipy
            maxclipy=tmpcl
        if minclipx > maxclipx:
            print('you switched min max in crop coordinates (longitude). fixing')
            tmpcl = minclipx
            minclipx=maxclipx
            maxclipx=tmpcl
        xrpha = xrpha.sel(x=slice(minclipx, maxclipx), y=slice(maxclipy, minclipy))
    if tolonlat:
        xrpha = xrpha.rename({'x': 'lon','y': 'lat'})
    return xrpha


def detrend_ifg_xr(xrda, isphase=True, return_correction = False, maxfringes = 4):
    """Estimates ramp of (wrapped) interferogram and corrects it. Based on Doris InSARMatlab Toolbox
    
    Args:
        xrda (xarray.Dataarray): input data array (interferogram)
        isphase (boolean): input array is phase (if not, expect complex ifg)
        return_correction (boolean): returns also the correction
        maxfringes (int): max amount of fringes to consider as proper correction
    
    Returns:
        xarray.Dataarray: dataarray of corrected ifg
    """
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

    
def filter_ifg_ml(ifg_ml, calc_coh_from_delta = False, radius = 1000, trunc = 4): #, sigma = 1, trunc = 2):  #, rotate = False):
    """Normalises interferogram and performs Gaussian filtering (expects proper structure of the ifg dataset).
    
    Args:
        ifg_ml (xarray.Dataset): input xr dataset (interferogram) - must contain ``pha`` data_var
        calc_coh_from_delta (boolean): will calculate local variance and use to improve ``gauss_coh`` measure
        radius (float): length of the Gaussian window in metres
        trunc (int): truncation of std dev for Gaussian window, by default trunc=4 and this is recommended for the shape
    Returns:
        xarray.Dataset: dataset that includes filtering results (as ``gauss_pha``, ``gauss_coh``, ``gauss_cpx``)
    """
    # get sigma, trunc from radius [m], converted to pixels using resolution
    resolution = get_resolution(ifg_ml, in_m=True)
    radius_px = radius/resolution
    #width_filter = 2*int( trunc*sigma + 0.5) +1  # definition within scipy ndimage filters.py - gaussian_filter1d
    #width_filter = 2*radius_px
    sigma = (radius_px - 1.5)/trunc
    #normalise mag
    tempar_mag1 = np.ones_like(ifg_ml.pha)
    ifg_ml['cpx'].values = magpha2RI_array(tempar_mag1, ifg_ml.pha.values)
    print('filter using (adapted) gauss filter')
    ifg_ml['gauss_cpx'] = filter_cpx_gauss(ifg_ml, sigma = sigma, trunc = trunc)
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



def wrap2phase(A):
    """Wraps array to -pi,pi (or 0,2pi?)
    """
    return np.angle(np.exp(1j*A))


def make_avg_amp(mlitiflist, hgtxr):
    """Generates average amplitude from list of MLI tiffs
    """
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
    """Generates standard deviation of amplitude from list of MLI tiffs
    """
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


def make_avg_coh(group, hgtxr):
    """Generates average coherence from 'group'
    """
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
    """Generates standard deviation of coherence from 'group'
    """
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
    """Builds amplitude stability map (or just avg/std amplitude) of a frame
    """
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
    """Builds coherence stability map (or just avg/std coherence) of a frame
    """
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






#### extra functions (not used in the workflow, but considered useful!)

def deramp_unw(xrda, dim=['lat','lon']):
    """Deramps unwrapped interferogram
    """
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
    """Deramps wrapped interferogram geotiff
    """
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
    
