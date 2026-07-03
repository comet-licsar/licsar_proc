#!/usr/bin/env python

import h5py
import datetime as dt
from datetime import datetime,timedelta
import numpy as np
import sys, os

#Usage:
#inh5 = 'S1C_SLC_066-0655-IW1-VV_20251123T052759.nc'
#create_TOPS_par_slab(inh5)


#try:
#    inh5 = sys.argv[1]
#    create_TOPS_par_slab(inh5)
#except:
#    print('Usage: gafa_licsar.py S1C_SLC_066-0655-IW1-VV_20251123T052759.nc')
#    exit()


'''
we need first to find raw files covering the frame coords (sent to MN),
then focus and merge into $epoch.IWx.slc files (and generate par files),
and then merge all to one frame SLC to be processed afterwards:
e.g. folders GAFA_1 GAFA_2
each has $epoch.IW?.slc* files inside

epoch=20250324
prevtab=''
slctab=$epoch.tab
mergetab=$epoch.merged.tab
createSLCtab $epoch merged.slc > $mergetab
createSLCtab $epoch slc > $slctab
for fld in GAFA_*; do
 createSLCtab $fld/$epoch slc > $fld.tab
 if [ ! -z $prevtab ]; then
   SLC_cat_ScanSAR $prevtab $fld.tab $mergetab
   for x in `ls $epoch.IW?.merged.slc*`; do mv $x `echo $x | sed 's/\.merged//'`; done
   prevtab=$slctab
 else
   prevtab=$fld.tab
 fi 
done
rm $mergetab
'''

'''
import LiCSquery as lq

lq.get_frame_files_date(frame, dt.datetime(2025,3,24), True

this will return list of SLC files registered within the frame (i.e. having frame bursts
'''


def identify_raw_covering_slc(slcfile, rawfiles):
    ''' Helper function to identify which of the raw files cover given slc file

    Args:
        slcfile: slc file identifier, e.g.
        rawfiles: list of raw filenames, e.g. [

    Returns:
        raw file id that covers the given slcfile
    '''
    t1slc, t2slc = get_dates_from_filename(slcfile)
    for rawfile in rawfiles:
        t1, t2 = get_dates_from_filename(rawfile)
        if t1 < t1slc and t2 > t2slc:
            return rawfile
    print('None of the raw files cover this SLC')
    return None


def get_dates_from_filename(filename):
    ''' Helper function to get t1 and t2 from given S1 filename

    Args:
        filename:  filename, e.g. 'S1A_IW_RAW__0SDV_20250324T233145_20250324T233217_058453_073AF4_530B'

    Returns:
        t1, t2:  dt.datetime
    '''
    t1, t2 = filename.split('_')[5:7]
    t1 = dt.datetime.strptime(t1, "%Y%m%dT%H%M%S")
    t2 = dt.datetime.strptime(t2, "%Y%m%dT%H%M%S")
    return t1, t2


def create_TOPS_par_slab(inh5):
    ''' Function to generate TOPS_par slab for given burst, e.g.

    Args:
        inh5: nc file for a burst, as output from GAFA. e.g. S1C_SLC_066-0655-IW1-VV_20251123T052759.nc

    '''
    f = h5py.File(inh5, "r")

    pos = f["/meta/grid/refsys/trajectory/pos"][:]
    vel = f["/meta/grid/refsys/trajectory/vel"][:]
    t   = f["/meta/grid/refsys/trajectory/utc"][:]

    units = f["/meta/grid/refsys/trajectory/utc"].attrs["units"].decode()
    ref = datetime.strptime(units.split("since")[1].strip(),
                            "%Y-%m-%d %H:%M:%S")
    times = [ref + timedelta(seconds=float(x)) for x in t]

    i = len(t)//2

    burst_date = times[0]
    burst_start_time = t[0]

    sensing_date = times[0]       # usually aligned with burst start
    sensing_start_time = t[0]

    print(f"burst_date_1:          {burst_date.isoformat()}")
    print(f"burst_start_time_1:    {burst_start_time:12.6f}  s")

    print(f"sensing_date_1:        {sensing_date.isoformat()}")
    print(f"sensing_start_time_1:  {sensing_start_time:12.6f}  s")

    print("sensing_position_1:", *pos[i])
    print("sensing_velocity_1:", *vel[i])

    dop_t = f["/meta/doppler/centroid/utc"][:]
    dop_units = f["/meta/doppler/centroid/utc"].attrs["units"].decode()

    ref = dop_units.split("since")[1].strip()
    t0 = datetime.strptime(ref.split(".")[0], "%Y-%m-%d %H:%M:%S")

    doppler_date = t0 + timedelta(seconds=float(dop_t[0]) * 1e-9)  # nanoseconds

    rgtime = f["/meta/doppler/centroid/rgtime"][:]
    doppler_srdelay = np.mean(rgtime)

    dop = f["/meta/doppler/centroid/values"][:]

    # Fit 4th-order polynomial
    coeffs = np.polyfit(rgtime, dop[0], 4)[::-1]

    # pad to 5 terms
    coeffs = np.pad(coeffs, (0, 5 - len(coeffs)))

    print(f"doppler_date_1:        {doppler_date.isoformat()}")
    print(f"doppler_time_1:        {doppler_srdelay:12.6f}   s")
    print(f"doppler_srdelay_1:     {doppler_srdelay:.5e}   s")

    print("doppler_polynomial_1:  " +
          " ".join(f"{c: .5e}" for c in coeffs))


    print("NOTE: The Azimuth FM is only approximated from velocity + geometry....")
    # crude approximation: use Doppler curvature
    az_fmrate_poly = coeffs * (-1e2)   # scaled derivative-like estimate

    print(f"az_fmrate_date_1:      {doppler_date.isoformat()}")
    print(f"az_fmrate_time_1:      {doppler_srdelay:12.6f}   s")
    print(f"az_fmrate_srdelay_1:   {doppler_srdelay:.5e}   s")

    print("az_fmrate_polynomial_1: " +
          " ".join(f"{c: .5e}" for c in az_fmrate_poly))

    grid_shape = f["/meta/grid/shape"][:]   # [lines, samples]
    gate_shape = f["/meta/gate/shape"][:]   # valid region

    lines, samples = grid_shape
    valid_lines, valid_samples = gate_shape


    first_valid_sample = 0
    last_valid_sample  = valid_samples

    first_valid_line   = 0
    last_valid_line    = valid_lines

    print(f"first_valid_sample_1:  {first_valid_sample:6d}         samples")
    print(f"last_valid_sample_1:   {last_valid_sample:6d}         samples")
    print(f"first_valid_line_1:    {first_valid_line:6d}         lines")
    print(f"last_valid_line_1:     {last_valid_line:6d}         lines")

    print(f"burst_boffset_1:       0")

    print("NOTE: safe approximation of ascending node fields")

    asc_node_delta = t[1] - t[0]   # time spacing
    burst_asc_node = 0.0           # placeholder

    print(f"asc_node_delta_1:      {asc_node_delta:12.6f}  s")
    print(f"burst_asc_node_1:      {burst_asc_node:12.6f}  bursts")


    print("NOTE: now only simple burst window estimates from timing and samples..")

    range_start = rgtime.min() * 1e6
    range_end   = rgtime.max() * 1e6

    az_start = t[0]
    az_end   = t[-1]

    print("burst_win_1:  "
          f"{range_start:.5f}   {range_end:.5f}   "
          f"{az_start:.6f}   {az_end:.6f}   "
          f"0   {samples}   1   {lines}")
    print("ext_burst_win_1:  "
          f"{range_start:.5f}   {range_end*1.001:.5f}   "
          f"{az_start:.6f}   {az_end*1.001:.6f}   "
          f"0   {samples}   1   {lines}")
    print("NEXT STEPS:")
    print(" - true TOPS burst segmentation")
    print(" - precise azimuth FM rates")
    print(" - exact azimuth timing grid")
    print(" - correct burst synchronization")
    f.close()


import glob
import xarray as xr

def merge_full(ncpath = os.getcwd()):
    ''' This will merge all nc burst files per each swath into merged.IWx.slc files, in FCOMPLEX format (2x float32)
    It now assumes all have the same number of samples, otherwise this will fail.
    '''
    curpth = os.getcwd()
    os.chdir(ncpath)
    ncs = glob.glob('S1?_SLC_*.nc')
    if not ncs:
        print('no slc files found')
        return False
    nc=ncs[0]
    datestr=nc.split('_')[-1].split('T')[0]
    for iw in ['IW1','IW2','IW3']:
        nclist = glob.glob('S1?_SLC_*-'+iw+'-*.nc')
        if nclist:
            nclist.sort()
            merge_ncs(nclist, outfile = datestr + '.'+iw+'.slc')
    os.chdir(curpth)


def _fix_shape_rg(cpx, dimrg):
    a,r = cpx.shape
    if r==dimrg:
        return cpx
    return np.pad(cpx, ((0, 0), (0, dimrg-r)),
                    mode='constant', constant_values = 0)


def merge_ncs(nclist, outfile = 'output.slc'):
    if os.path.exists(outfile):
        os.system('rm '+outfile)
    dimaz, dimrg = 0, 0
    for nc in nclist:
        a=xr.open_dataset(nc)
        dimaznc, dimrgnc = a.raster.shape
        #dimaz = max(dimaz,dimaznc)
        dimrg = max(dimrg,dimrgnc)
    print('Creating '+outfile+' in FCOMPLEX and dimensions of')
    print('samples: '+str(dimrg))
    # print('lines: '+str(int(dimaz*len(nclist))))
    # print('(this means '+str(dimrg)+' samples x '+str(dimaz)+' lines per burst)')
    for nc in nclist:
        print(nc)
        a=xr.open_dataset(nc)
        arr=a.raster.values
        cpx = (arr['r'] + 1j * arr['i']).astype(np.complex64)
        # cpx = _fix_shape(cpx, dimaz, dimrg)
        cpx = _fix_shape_rg(cpx, dimrg)
        cpx = cpx.byteswap()
        with open(outfile, 'ab') as f:
            cpx.tofile(f)


from pathlib import Path
import struct
import numpy as np


def load_gdr(fname):
    """ by Copilot
    Load a GDR raster and return. ifgs work ok, sigma is todo (not yet)

        array, metadata

    Supports complex64 and float32 rasters.

    e.g.
ifg, meta = load_gdr(
    "ifgm_2025-03-27_2025-04-02_143-0095-IW1-VV.gdr"
)
# MAG:
plt.imshow(np.log1p(np.abs(ifg)), cmap="gray")
plt.show()
# PHA:
plt.imshow(np.angle(ifg), cmap="hsv")
plt.colorbar()
plt.show()

    """
    fname = Path(fname)
    with open(fname, "rb") as f:
        # ---- file header ----
        magic = f.read(4)
        if magic[1:] != b"GDR":
            raise ValueError("Not a GDR file")
        version = struct.unpack("<I", f.read(4))[0]
        schema_len = struct.unpack("<Q", f.read(8))[0]
        # ---- schema text ----
        schema = f.read(schema_len).decode(
            "latin1",
            errors="ignore"
        )
        # infer raster dtype
        if "dtype:complex64" in schema:
            dtype = np.complex64
        elif "dtype:float32" in schema:
            dtype = np.float32
        elif "dtype:float64" in schema:
            dtype = np.float64
        else:
            raise ValueError("Unsupported raster dtype")
        # ---- beginning of metadata values ----
        meta_start = 16 + schema_len
        f.seek(meta_start)
        meta = f.read(512)
    # ------------------------------------------------------------------
    # Heuristic extraction of image shape
    #
    # We observed that the first grid.shape values appear early
    # within the metadata as two int64 values.
    # ------------------------------------------------------------------

    shape = None
    ints = np.frombuffer(meta, dtype="<i8")
    for i in range(len(ints) - 1):
        a = ints[i]
        b = ints[i + 1]
        if (
            1 < a < 100000
            and 1 < b < 100000
        ):
            pixels = a * b
            filesize = fname.stat().st_size
            data_bytes = pixels * np.dtype(dtype).itemsize
            if filesize > data_bytes:
                diff = filesize - data_bytes
                # metadata size should be plausible
                if 1000 < diff < 1000000:
                    shape = (int(a), int(b))
                    break
    if shape is None:
        raise RuntimeError("Could not determine raster shape")
    # ---- locate raster from end of file ----
    filesize = fname.stat().st_size
    nbytes = np.prod(shape) * np.dtype(dtype).itemsize
    data_offset = filesize - nbytes
    arr = np.fromfile(
        fname,
        dtype=dtype,
        offset=data_offset
    ).reshape(shape)
    metadata = {
        "version": version,
        "shape": shape,
        "dtype": dtype,
        "schema": schema,
        "data_offset": data_offset,
    }
    return arr, metadata

'''
# 1. merge the files (this will pad only in range)
cd 20250324; python3 -c "from gafa_licsar import *; merge_full()"
cd ../20250405; python3 -c "from gafa_licsar import *; merge_full()"
# 2. move the slcs and par files to some new directory
# 3. coregister them
createSLCtab 20250324 slc > 20250324.tab
createSLCtab 20250405 slc > 20250405.tab
createSLCtab 20250405 rslc > 20250405R.tab
# the previous commands i tried:
# SLC_mosaic_ScanSAR 20250324.tab 20250324.slc 20250324.slc.par 20 4 - - # to use existing burst windows
# SLC_mosaic_ScanSAR 20250324.tab 20250324.gamma.slc 20250324.gamma.slc.par 20 4 1 - # to have gamma recalc burst windows
# multi_look 20250324.slc 20250324.slc.par 20250324.mli 20250324.mli.par 20 4
# but the idea is to coregister per swath... and then use existing routines to create bovl. so:
ScanSAR_coreg.py 20250324.tab 20250324 20250405.tab 20250405 20250405R_tab - 20 4 --r_only --no_cleaning # the following params might be useful to check/test: --no_check or --burst_check
# 4. check the interferograms in the temp folder, or create the bovls
'''