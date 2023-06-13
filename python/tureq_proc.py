#!/usr/bin/env python3

# LOADING FUNCTIONS
from lics_unwrap import *
import dask.array as da
#import dask_ndfilters as ndfilters    # use this if the below doesn't work!
from dask_image import ndfilters

from scipy.ndimage import generic_filter
#from skimage.morphology import disk

def get_disk_ones(block):
    mask=np.zeros(block.shape) #should be square
    #nyquistlen=int(mask.shape[0]/2+0.5) + 1 #+ extrapx
    circle=unit_circle(int(mask.shape[0]/2+0.5)-1) #will contain +1 px for zero
    i=int((mask.shape[0]-circle.shape[0])/2+0.5)
    j=int((mask.shape[1]-circle.shape[1])/2+0.5)
    mask[i:i+circle.shape[0],j:j+circle.shape[1]]=circle
    mask[mask==0]=np.nan
    return mask


def filterhistmed(block, amin, amax, bins=20, medbin=True):
    """Support function to be used with generic_filter (where only 1D array is passed, expecting one output->using median here only
    """
    if np.isnan(block).all():
        return np.nan
    histc, histe = np.histogram(block,range=(amin,amax), bins=bins)
    histmax=np.argmax(histc)
    #minval=histe[histmax]
    #maxval=histe[histmax+1]
    # tiny update
    if histmax == 0:
        histmax=1
    elif histmax == bins-1:
        histmax = bins-2
    minval=histe[histmax-1]
    maxval=histe[histmax+2]
    # add median or interpolate:
    if medbin:
        bb=block[block>minval]
        bb=bb[bb<maxval]
        outval = np.nanmedian(bb)
    else:
        try:
            blockxr=xr.DataArray(block)
            blockxr.values=interpolate_nans(blockxr.where(blockxr>minval).where(blockxr<maxval).values, method='linear')
            blockxr.values=interpolate_nans(blockxr.values, method='nearest') # just to be sure..
            outval = float(blockxr.values[int(block.shape[0]/2),int(block.shape[1]/2)])
        except:
            outval = np.nan
    return outval


def filter_histmed_ndarray(ndarr, winsize=32, bins=20, medbin=True):
    """Main filtering function, works with both numpy.ndarray and xr.DataArray
    Args:
        medbin (boolean): if False, it will interpolate (fit) the central value from the bin subset. otherwise returns its median
    """
    #footprint=disk(winsize)
    amin=np.nanmin(ndarr)
    amax=np.nanmax(ndarr)
    footprint=unit_circle(int(winsize/2))
    return generic_filter(ndarr, filterhistmed, footprint=footprint, mode='constant', cval=np.nan,
                      extra_keywords= {'amin': amin, 'amax':amax, 'bins':bins, 'medbin':medbin})

'''
from skimage.restoration import inpaint
def inpaintxr(xra):
    xrb=xra.copy()
    mask = np.isnan(xra.values)*1
    xrb.values = inpaint.inpaint_biharmonic(xra.values, mask)
    return xrb
'''

origfigsize=plt.rcParams['figure.figsize']
def plot2(A,B = None):
    if type(B) == type(None):
        numpl = 2
    else:
        numpl = 3
    plt.rcParams["figure.figsize"] = [int(6*numpl),4]
    plt.subplot(1,numpl,1)
    A.rename('px').plot()
    plt.subplot(1,numpl,2)
    A.plot.hist(bins=20)
    plt.axvline(A.median().values, color='black')
    #B.rename('rad').plot()
    if numpl > 2:
        plt.subplot(1,numpl,3)
        #AA.toremove
        B.plot()
    plt.show()
    plt.rcParams['figure.figsize']=origfigsize
    return


# older filter, working but.. median..
def medianfilter_array(arr, ws = 32):
    """use dask median filter on array
    works with both xarray and numpy array
    """
    chunksize = (ws*8, ws*8)
    if type(arr)==type(xr.DataArray()):
        inn = arr.values
    else:
        inn = arr
    arrb = da.from_array(inn, chunks=chunksize)
    arrfilt=ndfilters.median_filter(arrb, size=(ws,ws), mode='reflect').compute()
    if type(arr)==type(xr.DataArray()):
        out = arr.copy()
        out.values = arrfilt
    else:
        out = arrfilt
    return out

# very gross initial filter
def init_filter(prevest, ml=20):
    tak=prevest.copy()
    prevest = prevest.coarsen({prevest.dims[0]:ml,prevest.dims[1]:ml}, boundary='trim').median()
    prevest = medianfilter_array(prevest, ws = 8)
    #prevest.values = interpolate_nans(prevest,'linear')
    prevest.values = interpolate_nans(prevest,'nearest') # for pixels outside of interpolation
    prevest = medianfilter_array(prevest, ws = 8)
    #xr.DataArray(azinn).plot()
    #prevest = prevest.interp_like(tak, method='linear')
    prevest = prevest.interp_like(tak, method='nearest')
    prevest.values = interpolate_nans(prevest,'nearest')
    prevest = medianfilter_array(prevest, ws = 16)
    return prevest


def filterhist_value(block, amin, amax, bins=20, medbin=True):
    #histc, histe = np.histogram(block,range=(amin,amax), bins=bins)
    #histmax=np.argmax(histc)
    #minval=histe[histmax]
    #maxval=histe[histmax+1]
    #
    #if medbin:
    #    outval = float(blockxr.where(blockxr>minval).where(blockxr<maxval).median())
    #else:
    #    blockxr.values=interpolate_nans(blockxr.where(blockxr>minval).where(blockxr<maxval).values, method='nearest')
    #    outval = float(blockxr.sel(a=blockxr.shape[0]/2,r=blockxr.shape[1]/2,method='nearest'))
    #
    if np.isnan(block).all():
        return np.nan
    try:
        # dask may change window size to very small ones. for this, we will skip the disk
        block=get_disk_ones(block)*block
    except:
        pass
    histc, histe = np.histogram(block,range=(amin,amax), bins=bins)
    histmax=np.argmax(histc)
    minval=histe[histmax]
    maxval=histe[histmax+1]
    #
    if medbin:
        bb=block[block>minval]
        bb=bb[bb<maxval]
        outval = np.nanmedian(bb)
    else:
        blockxr=xr.DataArray(block)
        try:
            blockxr.values=interpolate_nans(blockxr.where(blockxr>minval).where(blockxr<maxval).values, method='linear')
            blockxr.values=interpolate_nans(blockxr.values, method='nearest') # just to be sure..
            outval = float(blockxr.values[int(block.shape[0]/2),int(block.shape[1]/2)])
        except:
            outval = np.nan
    return outval


import time

_start_time = time.time()


def tic():
    global _start_time
    _start_time = time.time()


def tac():
    t_sec = round(time.time() - _start_time)
    (t_min, t_sec) = divmod(t_sec, 60)
    (t_hour, t_min) = divmod(t_min, 60)
    print('Elapsed time: {}h:{}m:{}s'.format(t_hour, t_min, t_sec))


def calculate_gradient(xar, deramp=False):
    """Calculates gradient of continuous data (not tested for phase)

    Args:
        xar (xr.DataArray): e.g. ifg['unw']
        deramp (bool): if True, it will remove overall ramp

    Returns:
        xr.DataArray
    """
    gradis = xar.copy()
    vgrad = np.gradient(gradis.values)
    gradis.values = np.sqrt(vgrad[0] ** 2 + vgrad[1] ** 2)
    if deramp:
        gradis = deramp_unw(gradis)
    return gradis


# use:
#    tic()
#    f()
#    tac()

















# LOADING DATA
#subdir='/gws/nopw/j04/nceo_geohazards_vol1/projects/LiCS/proc/current/subsets/test_tur_rs/021D'
subdir='/work/scratch-pw3/licsar/earmla/batchdir/021D_05266_252525'
offs=os.path.join(subdir,'OFF/20230129_20230210/tracking.offsets') # or just remove the '128x128'
# outputs:
outcpxfile=os.path.join(subdir,'OFF/20230129_20230210/tracking.offsets.filtered')
outlutfile=os.path.join(subdir,'OFF/20230129_20230210/offsets.filtered.lut')
outlutfilefull=os.path.join(subdir,'OFF/20230129_20230210/offsets.filtered.lut.full')

#lena=1409
#lenr=1278
# see OFF/20230129_20230210/tracking.off:
lena=8662
lenr=3447
#outlenrng=25561
#outlenazi=5636
# see RSLC/20230129/20230129.rslc.par
outlenrng=68953
outlenazi=34651
thresm=5.5
rngres=2.329562 #precision will be needed for range
azires=14

offs = np.fromfile(offs, dtype=np.complex64).byteswap().reshape((lena,lenr))
#rng = os.path.join(subdir,'IFG/20230129_20230210/disp_map.rng')
#azi = os.path.join(subdir,'IFG/20230129_20230210/disp_map.azi')
#corr=os.path.join(subdir,'OFF/20230129_20230210.128x128/tracking.corr')
#rng = np.fromfile(rng, dtype=np.float32).byteswap().reshape((1409,1278))
#azi = np.fromfile(azi, dtype=np.float32).byteswap().reshape((1409,1278))
#corr = np.fromfile(corr, dtype=np.float32).byteswap().reshape((lena,lenr))
# remove values above thresm
#thresm=5
#rng = rng[np.abs(rng)<thresm]
#azi = azi[np.abs(azi)<thresm]

rng=xr.DataArray(np.real(offs))
azi=xr.DataArray(np.imag(offs))
#rng = rng.expand_dims({"x": rng.dim_1, "y": rng.dim_0})
rng = xr.DataArray(
    data=rng.values,
    dims=["a", "r"],
    coords={"a": rng.dim_0.values, "r": rng.dim_1.values},
)
azi = xr.DataArray(
    data=azi.values,
    dims=["a", "r"],
    coords={"a": azi.dim_0.values, "r": azi.dim_1.values},
)
#azi=azi.expand_dims({"x": azi.dim_1, "y": azi.dim_0})
#corr=xr.DataArray(corr)
#corr = xr.DataArray(
#    data=corr.values,
#    dims=["a", "r"],
#    coords={"a": corr.dim_0.values, "r": corr.dim_1.values},
#)
#corr=corr.expand_dims({"x": corr.dim_1, "y": corr.dim_0})
rng=rng.where(rng!=0)
azi=azi.where(azi!=0)

rng=rng.where(np.abs(rng)<thresm/rngres)
azi=azi.where(np.abs(azi)<thresm/azires)
# but also those cross-errors:
rng=rng.where(np.abs(azi)<thresm/azires)
azi=azi.where(np.abs(rng)<thresm/rngres)

rng.values=remove_islands(rng.values, pixelsno = 25)
azi.values=remove_islands(azi.values, pixelsno = 25)


print('filtering azi')
tic()
azifilt = filter_histmed_ndarray(azi, winsize=128, bins=10)
medres2 = (azi-azifilt).copy()
medres2=medres2.fillna(0)
medres2=medianfilter_array(medres2, ws=64)
tac()

print('filtering rng')
tic()
rngfilt = filter_histmed_ndarray(rng, winsize=64, bins=10)
medres = (rng-rngfilt).copy()
medres=medres.fillna(0)
medres=medianfilter_array(medres, ws=32)
outrng=rngfilt+medres
tac()

# STORING DATA

fullazi=outazi
fullrng=outrng

# store back to cpx
#outcpxfile='outcpxfile'
outcpx = fullrng.values + 1j* fullazi.values
outcpx.astype('complex64').tofile(outcpxfile)

# but we need LUT, so store by:
# a) resampling to whole dimensions:
z = xr.DataArray(np.zeros(outlenrng*outlenazi).reshape(outlenazi, outlenrng),
                 dims=["a", "r"],
                 coords={"a": np.arange(outlenazi), "r": np.arange(outlenrng)}
                )

# first, export the LUT as is (multilooked)
if np.max(np.isnan(fullrng)):
    method = 'nearest'
    fullrng.values = interpolate_nans(fullrng.values, method=method)


if np.max(np.isnan(fullazi)):
    method = 'nearest'
    fullazi.values = interpolate_nans(fullazi.values, method=method)


# adding the pixel numbers themselves:
rnglut = fullrng + np.tile(fullrng.r.values, (len(fullrng.a.values),1))
azilut = fullazi + np.tile(fullazi.a.values, (len(fullazi.r.values),1)).T
# storing
print('storing')
outlut = rnglut.values + 1j* azilut.values
#outlutfile='outlutfile'
outlut.astype('complex64').byteswap().tofile(outlutfile)

# store also the range offsets to be used as prevest
rnginmm = fullrng.values*rngres*1000
mm2rad_s1(rnginmm).astype('float32').tofile('rngoffsets_prevest_LE')


# a.2) full res output
# for interpolation, i need to have the corresponding px set. hope i did not miss +1 pixel...
fullrngfull=fullrng.assign_coords(a=fullrng.a.values*int(outlenazi/lena), r=fullrng.r.values*int(outlenrng/lenr))
fullazifull=fullazi.assign_coords(a=fullazi.a.values*int(outlenazi/lena), r=fullazi.r.values*int(outlenrng/lenr))
print('full interpolation to original RSLC dimensions')
fullrngfull = fullrngfull.interp_like(z,method='nearest',kwargs={"fill_value": "extrapolate"})
fullazifull = fullazifull.interp_like(z,method='nearest',kwargs={"fill_value": "extrapolate"})
if np.max(np.isnan(fullrngfull)):
    method = 'nearest'
    fullrngfull.values = interpolate_nans(fullrngfull.values, method=method)


if np.max(np.isnan(fullazifull)):
    method = 'nearest'
    fullazifull.values = interpolate_nans(fullazifull.values, method=method)


# b) adding the pixel numbers themselves:
rnglut = fullrngfull + np.tile(fullrngfull.r.values, (len(fullrngfull.a.values),1))
azilut = fullazifull + np.tile(fullazifull.a.values, (len(fullazifull.r.values),1)).T
# c) storing
print('storing')
outlut = rnglut.values + 1j* azilut.values
#outlutfile='outlutfile'
outlut.astype('complex64').byteswap().tofile(outlutfilefull)

print('done')




# and finally, recreate the RSLC with that:
#cd /gws/nopw/j04/nceo_geohazards_vol1/projects/LiCS/proc/current/subsets/test_tur_rs/021D
#mv ~/outlutfile OFF/20230129_20230210/offsets.filtered.lut.full
# now I should be able to do in GAMMA:
# SLC_interp_lt RSLC/20230210.orig/20230210.rslc RSLC/20230129.orig/20230129.rslc.par RSLC/20230129.orig/20230129.rslc.par OFF/20230129_20230210/offsets.filtered.lut.full RSLC/20230129.orig/20230129.rslc.par RSLC/20230210.orig/20230210.rslc.par - RSLCRS/20230210/20230210.rslc RSLCRS/20230210/20230210.rslc.par - - 5
# or, using :
mkdir -p temp RSLCRS/$s
SLC_interp_lt_ScanSAR tab/$s'R_tab' RSLC/$s/$s.rslc.par tab/$m'R_tab' RSLC/$m/$m.rslc.par OFF/20230129_20230210/offsets.filtered.lut.full RSLC/$m/$m.rslc.par RSLC/$s/$s.rslc.par - tab/$s'_RS RSLCRS/$s/$s.rslc RSLCRS/$s/$s.rslc.par - 5 temp

# this way i resample the 20230210 using the fine fullres offsets in the LUT to the new RSLC
cd RSLCRS/20230210; multi_look 20230210.rslc 20230210.rslc.par 20230210.rslc.mli 20230210.rslc.mli.par 20 4; cd ../..
# echo $m'_'$
# rm -r IFG/20230129_20230210; LiCSAR_03_mk_ifgs.py -d . -i ifg.list -a 4 -r 20
