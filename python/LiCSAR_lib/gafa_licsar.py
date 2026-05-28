#!/usr/bin/env python

import h5py
import datetime as dt
from datetime import datetime,timedelta
import numpy as np
import sys

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
for fld in GAFA_1 GAFA_2; do
 createSLCtab GAFA_1/$epoch slc > $fld.tab
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
    for iw in ['IW1','IW2','IW3']:
        nclist = glob.glob('S1?_SLC_*-'+iw+'-*.nc')
        if nclist:
            nclist.sort()
            merge_ncs(nclist, outfile = 'merged.'+iw+'.slc')
    os.chdir(curpth)


def merge_ncs(nclist, outfile = 'output.slc'):
    if os.path.exists(outfile):
        os.system('rm '+outfile)
    for nc in nclist:
        print(nc)
        a=xr.open_dataset(nc)
        arr=a.raster.values
        cpx = (arr['r'] + 1j * arr['i']).astype(np.complex64)
        cpx = cpx.byteswap()
        with open(outfile, 'ab') as f:
            cpx.tofile(f)
