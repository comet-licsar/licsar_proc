#!/usr/bin/env python3

# Milan Lazecky, 2023
# not-optimized approach for rubber-sheeting, as applied to the coseismic data of 2023 Turkey-Syria earthquake (not published yet, planned to SARWatch)
# contains specific filtering that might be improved but is better edge preserving than other simple approaches, such as median filter


print('''
#  How I created tracking offsets? after few iterations, i end up in this way (note deramp is ON here, for TOPS, and note multilooking factor 20/4!):
m=20230129
s=20230210
outdir=OFF/$m'_'$s
mkdir -p $outdir
create_offset RSLC/$m/$m.rslc.par RSLC/$s/$s.rslc.par $outdir/tracking.off 1 20 4 0
time offset_pwr_tracking RSLC/$m/$m.rslc RSLC/$s/$s.rslc RSLC/$m/$m.rslc.par RSLC/$s/$s.rslc.par $outdir/tracking.off $outdir/tracking.offsets $outdir/tracking.corr 128 128 - 2 0.1 20 4 - - - - - - 1 - - - $outdir/ccs
# result should be SLC2 offsets towards SLC1 (S->M)
''')

# IMPORTS
from lics_unwrap import *
import dask.array as da
from dask_image import ndfilters
from scipy.ndimage import generic_filter
import time
import os

# SETTINGS
full = False # for testing of 1:1 LUT table. Do not use/not tested/not working
filt_hist = True # histogram-based edge-preserving filter. SUPER SLOW but much better than median-based (will run if this is False)

subdir= os.getcwd() #'/work/scratch-pw3/licsar/earmla/batchdir/021D_05266_252525'
m='20230129'
s='20230210'
lena=1409
lenr=1278
mlrng=20  # this should be same as rstep in offset_pwr_tracking
mlazi=4
# see OFF/20230129_20230210/tracking.off
if full:
    outlenrng=68953
    outlenazi=34651
    # see RSLC/20230129/20230129.rslc.par

thresm=5.5  # will cut-off offsets larger than this, in metres
rngres=2.329562 #precision will be needed for range
azires=14


####################### Nothing more to set, just run this and see the message below

# LOADING FUNCTIONS

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






# LOADING/PROCESSING PART

offdir = os.path.join(subdir,'OFF',m+'_'+s)
offs=os.path.join(offdir,'tracking.offsets') # or just remove the '128x128'
# outputs:
outcpxfile=os.path.join(offdir,'tracking.offsets.filtered')
outlutfile=os.path.join(offdir,'offsets.filtered.lut')

if full:
    outlutfilefull=os.path.join(offdir,'offsets.filtered.lut.full')


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

# check the output as e.g.
# rng.plot(vmax=1); plt.show()

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
if filt_hist:
    azifilt = filter_histmed_ndarray(azi, winsize=128, bins=10)
else:
    azifilt=medianfilter_array(azi, ws=64)


medres = (azi-azifilt).copy()
medres=medres.fillna(0)
medres=medianfilter_array(medres, ws=64)
outazi=azifilt+medres
tac()

print('filtering rng')
tic()
if filt_hist:
    rngfilt = filter_histmed_ndarray(rng, winsize=64, bins=10)
else:
    rngfilt=medianfilter_array(rng, ws=32)


medres = (rng-rngfilt).copy()
medres=medres.fillna(0)
medres=medianfilter_array(medres, ws=32)
outrng=rngfilt+medres
tac()

# STORING DATA

# store back to cpx
#outcpxfile='outcpxfile'
outcpx = outrng.fillna(0).values + 1j* outazi.fillna(0).values
outcpx.astype('complex64').tofile(outcpxfile)

# but we need LUT, so:
# first, export the LUT as is (multilooked)
if np.max(np.isnan(outrng)):
    print('warning, nans remain - interpolating by NN')
    method = 'nearest'
    outrng.values = interpolate_nans(outrng.values, method=method)


if np.max(np.isnan(outazi)):
    print('warning, nans remain - interpolating by NN')
    method = 'nearest'
    outazi.values = interpolate_nans(outazi.values, method=method)

# 2023/06/30 - update!! we can actually use SLC_interp_map !!
#outcpx = outrng.fillna(0).values + 1j* outazi.fillna(0).values
#outcpx.astype('complex64').byteswap.tofile(outcpxfile+'.BE4')


# need to fix the lut for multilook factor!
# adding the pixel numbers themselves:
rnglut = outrng/mlrng + np.tile(outrng.r.values, (len(outrng.a.values),1))
azilut = outazi/mlazi + np.tile(outazi.a.values, (len(outazi.r.values),1)).T

# storing
print('storing')
outlut = rnglut.values + 1j* azilut.values
#outlutfile='outlutfile'
outlut.astype('complex64').byteswap().tofile(outlutfile)

# store also the range offsets to be used as prevest - see lics_unwrap.py
rnginmm = outrng.values*rngres*1000
mm2rad_s1(rnginmm).astype('float32').tofile('rngoffsets_prevest_LE')


if full:
    # a) resampling to whole dimensions:
    z = xr.DataArray(np.zeros(outlenrng*outlenazi).reshape(outlenazi, outlenrng),
                 dims=["a", "r"],
                 coords={"a": np.arange(outlenazi), "r": np.arange(outlenrng)}
                )
    # a.2) full res output
    # for interpolation, i need to have the corresponding px set. hope i did not miss +1 pixel...
    fullrngfull=outrng.assign_coords(a=outrng.a.values*int(outlenazi/lena), r=outrng.r.values*int(outlenrng/lenr))
    fullazifull=outazi.assign_coords(a=outazi.a.values*int(outlenazi/lena), r=outazi.r.values*int(outlenrng/lenr))
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


print('done, see below for using the offsets-based LUT that was stored as:')
print(outlutfile)

print('''
finally, resample the RSLC using the updated offsets-based LUT file:
m=20230129
s=20230210
offlut=OFF/$m'_'$s/offsets.filtered.lut
slc2=RSLC/$s/$s.rslc
slc1=RSLC/$m/$m.rslc
outdir=RSLCRS
mkdir -p $outdir/$s
# indeed, SLC_interp_lt resamples slc2->slc1, so the offsets should correspond to this
SLC_interp_lt $slc2 $slc1.par $slc2.par $offlut $slc1.mli.par $slc2.mli.par - $outdir/$s/$s.rslc $outdir/$s/$s.rslc.par - - 5
# or:
# offpar=OFF/$m'_'$s/tracking.off
# swap_bytes OFF/$m'_'$s/tracking.offsets.filtered OFF/$m'_'$s/tracking.offsets.filtered.BE4 4
# interp_ad OFF/$m'_'$s/tracking.offsets.filtered.BE4  OFF/$m'_'$s/tracking.offsets.filtered.BE4.interpolated 1278 - - - - 0
# SLC_interp_map $slc2 $slc1.par $slc2.par $offpar $outdir/$s/$s.rslc $outdir/$s/$s.rslc.par OFF/$m'_'$s/tracking.off OFF/$m'_'$s/tracking.offsets.filtered.BE4.interpolated

# and.... IT WORKED!!!!!
# and with the LUT it is also OK!!!! there is a tiny diff (not sure why) but the coh was improved in almost same way (0+-0.01)

# so for the burst ovls, i can do following (checked, works, not sure if any real effect):
# SLC_interp_lt_ScanSAR tab/20230210R_tab RSLC/20230210/20230210.rslc.par tab/20230129_tab RSLC/20230129/20230129.rslc.par OFF/20230129_20230210/offsets.filtered.lut RSLC/20230129/20230129.rslc.mli.par RSLC/20230210/20230210.rslc.mli.par OFF/20230129_20230210/tracking.off tab/20230210rsR_tab rsRSLC/20230210/20230210.rslc rsRSLC/20230210/20230210.rslc.par


cd $outdir/$s; multi_look $s.rslc $s.rslc.par $s.rslc.mli $s.rslc.mli.par 20 4; cd ../..

# now it is possible to generate interferograms and check, e.g.
echo $m'_'$s > ifg.list
rm -r IFG/$m'_'$s 2>/dev/null; mkdir -p IFG
if [ ! -d RSLC/$s.orig ]; then
 cd RSLC; mv $s $s.orig; ln -s ../$outdir/$s; cd ..
fi
mkdir -p log
LiCSAR_03_mk_ifgs.py -d . -i ifg.list -a 4 -r 20
rm RSLC/$s; mv RSLC/$s.orig RSLC/$s
display IFG/20230129_20230210/20230129_20230210.diff.bmp &


# i tested also max Lanczos degree: 9. No effect for coherence
''')








'''

some older notes here, just in case if useful for the full approach (that failed):

# and finally, recreate the RSLC with that - first in the subset:
cd /gws/nopw/j04/nceo_geohazards_vol1/projects/LiCS/proc/current/subsets/test_tur_rs/021D
#mv ~/outlutfile OFF/20230129_20230210/offsets.filtered.lut
#mv ~/outlutfilefull OFF/20230129_20230210/offsets.filtered.lut.full

# with the small LUT:
slc2=RSLC/20230210/20230210.rslc
slc1=RSLC/20230129/20230129.rslc
mli2=RSLC/20230210/20230210.rslc.mli
SLC_interp_lt $slc2 $slc1.par $slc2.par OFF/20230129_20230210/offsets.filtered.lut $slc1.mli.par $slc2.mli.par - RSLCRS2/20230210/20230210.rslc RSLCRS2/20230210/20230210.rslc.par - - 5
cd RSLCRS2/20230210; multi_look 20230210.rslc 20230210.rslc.par 20230210.rslc.mli 20230210.rslc.mli.par 20 4; cd ../..

# or (better?):
mkdir tempp;
SLC_interp_lt_ScanSAR tab/20230210R_tab RSLC/20230210/20230210.rslc.par tab/20230129R_tab RSLC/20230129/20230129.rslc.par OFF/20230129_20230210/offsets.filtered.lut.small RSLC/20230129/20230129.rslc.mli.par RSLC/20230210/20230210.rslc.mli.par - tab/20230210RS_tab RSLCRS/20230210/20230210.rslc RSLCRS/20230210/20230210.rslc.par - 6 tempp


# to make ifgs
rm RSLC/20230210
ln -s `pwd`/RSLCRS2/20230210 `pwd`/RSLC/20230210
#SLC_interp_lt RSLC/20230210.orig/20230210.rslc RSLC/20230210.orig/20230210.rslc.par RSLC/20230210.orig/20230210.rslc.par OFF/20230129_20230210/offsets.filtered.lut RSLC/20230210.orig/20230210.rslc.mli.par RSLC/20230210.orig/20230210.rslc.mli.par - RSLC2/20230210/20230210.rslc RSLC2/20230210/20230210.rslc.par - - 8
#multi_look 20230210.rslc 20230210.rslc.par 20230210.rslc.mli 20230210.rslc.mli.par 20 4


# now I should be able to do in GAMMA:
# SLC_interp_lt RSLC/20230210.orig/20230210.rslc RSLC/20230129.orig/20230129.rslc.par RSLC/20230129.orig/20230129.rslc.par OFF/20230129_20230210/offsets.filtered.lut.full RSLC/20230129.orig/20230129.rslc.par RSLC/20230210.orig/20230210.rslc.par - RSLCRS/20230210/20230210.rslc RSLCRS/20230210/20230210.rslc.par - - 5
# or, using :
mkdir -p temp RSLCRS/$s
SLC_interp_lt_ScanSAR tab/$s'R_tab' RSLC/$s/$s.rslc.par tab/$m'R_tab' RSLC/$m/$m.rslc.par OFF/20230129_20230210/offsets.filtered.lut.full RSLC/$m/$m.rslc.par RSLC/$s/$s.rslc.par - tab/$s'_RS RSLCRS/$s/$s.rslc RSLCRS/$s/$s.rslc.par - 5 temp
# this is to use the not-full version
SLC_interp_lt_ScanSAR tab/$s'R_tab' RSLC/$s/$s.rslc.par tab/$m'R_tab' RSLC/$m/$m.rslc.par OFF/20230129_20230210/offsets.filtered.lut.full RSLC/$m/$m.rslc.par RSLC/$s/$s.rslc.par - tab/$s'_RS RSLCRS/$s/$s.rslc RSLCRS/$s/$s.rslc.par - 5 temppp

# this way i resample the 20230210 using the fine fullres offsets in the LUT to the new RSLC
cd RSLCRS/20230210; multi_look 20230210.rslc 20230210.rslc.par 20230210.rslc.mli 20230210.rslc.mli.par 20 4; cd ../..
# echo $m'_'$
# rm -r IFG/20230129_20230210; LiCSAR_03_mk_ifgs.py -d . -i ifg.list -a 4 -r 20
'''
