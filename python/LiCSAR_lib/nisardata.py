#!/usr/bin/env python

'''
Created by Milan Lazecky in 2026 - based on his previous dask experience.
Yes, I first tried Copilot and I keep his messy codes in captions here, for now..
Then I looked at my lics_unwrap, reminded few tricks, and replaced the many lines by few.. and more correct ones.
This is only to remind that Copilot is really useful for a fast first step, but you just need the experience and some digging.

This is still early in dev, but NISAR search, download, and ifg processing works well - give it a try!
(just get your GSLC .h5 files of the same track and path, and then apply generate_ifg function)
'''
import os
import datetime as dt
import re
import pandas as pd
import geopandas as gpd
import requests
from LiCSAR_misc import *
import asf_search as asf
import time
import lics_processing as lp
from shapely.geometry import shape
from shapely import wkt
from shapely.ops import transform
from pyproj import Transformer
import h5py
import dask.array as da
import xarray as xr
import numpy as np
from rasterio.crs import CRS  # optional but convenient

from dask.diagnostics import ProgressBar
ProgressBar().register()


def wgs2utm(polygon, target_crs):
    ''' the bbox should be shapely polygon in WGS-84'''
    source_crs = "EPSG:4326"
    transformer = Transformer.from_crs(source_crs, target_crs, always_xy=True)
    poly_utm = transform(transformer.transform, bbox)
    return poly_utm


def fullchain(lon1, lat1, lon2, lat2, nisarslcpath = '/gws/ssde/j25a/nceo_geohazards/vol1/public/shared/NISAR/allinputs', downloadit = True,
              clipit = True, processit = True):
    # will get automatically
    if not nisarslcpath:
        if 'XFCPATH' in os.environ:
            nisarslcpath = os.path.join(os.environ['XFCPATH'], 'SLC')
        elif 'LiCSAR_SLC' in os.environ:
            nisarslcpath = os.environ['LiCSAR_SLC']
        else:
            nisarslcpath = os.getcwd()
        print('expected input data into '+nisarslcpath)
    #
    bbox = lp.bbox_to_wkt(lon1, lat1, lon2, lat2)
    nsrs = get_nisar_data(bbox, outAspd = True)
    # filter a bit:
    polygon = wkt.loads(bbox)
    nsrs=nsrs[nsrs.intersects(polygon)]
    # keep only data needed for whole bbox:
    nsrsel=nsrs.groupby(['flightDirection','pathNumber', 'frameNumber']).head(1)[['flightDirection', 'pathNumber', 'frameNumber','geometry']]
    for i, ln in nsrsel[nsrsel.contains(polygon)].iterrows():
        todrop=nsrs[nsrs['flightDirection']==ln['flightDirection']][nsrs['pathNumber']==ln['pathNumber']][nsrs['frameNumber']!=ln['frameNumber']]
        nsrs=nsrs.drop(todrop.index)
    # update the selection
    nsrsel = nsrs.groupby(['flightDirection', 'pathNumber', 'frameNumber']).head(1)[['flightDirection', 'pathNumber', 'frameNumber', 'geometry']]
    # now can download them
    print('data available from '+str(len(nsrsel.groupby(['flightDirection','pathNumber'])))+' orbital passes')
    if downloadit:
        print('Now we check and download ' + str(len(nsrs)) + ' files')
        for i, ln in nsrs.iterrows():
            print('downloading '+ln['sceneName']+'.h5')
            url = ln['url']
            fpath = download(url, nisarslcpath)  # using this way, because for some reason i cannot get ASF session established..
            if not os.path.exists(fpath):
                print('some error in '+fpath)
    if processit:
        freq_code = 'B'
        if clipit:
            # need to reproject the bbox to given UTM.. will be done in load function
            clipping_bbox = bbox
        else:
            clipping_bbox = None
        nsrs=nsrs.sort_values('startTime')  # sort since the beginning
        for i, sset in nsrsel.iterrows():
            opass=sset['flightDirection']
            pan = sset['pathNumber']
            frn = sset['frameNumber']
            frame = opass[0]+'.'+str(pan)+'.'+str(frn)
            framedir = frame
            if not os.path.exists(framedir):
                os.mkdir(framedir)
            print('processing frame '+frame)
            tmpsel = nsrs[nsrs['flightDirection'] == opass][nsrs['pathNumber'] == pan][nsrs['frameNumber'] == frn]
            # now lets get network...
            ifgs = get_network(tmpsel, ntype='daisy')
            for ifg in ifgs:
                in1 = os.path.join(nisarslcpath, ifg[0].sceneName+'.h5')
                in2 = os.path.join(nisarslcpath, ifg[1].sceneName + '.h5')
                epoch1 = ifg[0].startTime.split('T')[0].replace('-','')
                epoch2 = ifg[1].startTime.split('T')[0].replace('-', '')
                if os.path.exists(in1) and os.path.exists(in2):
                    outncfile = os.path.join(framedir, epoch1+'_'+epoch2+'.nc')
                    generate_ifg(
                        in1=in1,
                        in2=in2,
                        freq_code=freq_code, polarization='HH', clipping_bbox=bbox,
                        target_resolution_m=110, outncfile=outncfile, create_wgs84_previews=True)
                else:
                    print('ERROR, file '+in1+' does not exist')


def get_network(tmpsel, ntype='daisy'):
    ifgs = []
    if ntype == 'daisy':
        for i in range(len(tmpsel)-1):
            ifgs.append([tmpsel.iloc[i], tmpsel.iloc[i+1]])
    else:
        print('only daisy chain implemented for now')
        return False
    return ifgs


# say we want to get NISAR data covering particular location, or region:
def get_nisar_data(wkt, dtype = 'GSLC', startdate = dt.datetime.strptime('20250101','%Y%m%d').date(),
             enddate = dt.date.today(), outAspd = False):
    ''' main search engine for NISAR
    you need to provide wkt as input - you can use lp.cliparea_geo2coords for this, or use e.g.
    lat=33.74; lon=-118.37
    wkt = f"POINT({lon} {lat})"
    '''
    #
    results=asf.geo_search(dataset='NISAR', processingLevel=dtype, intersectsWith=wkt, 
                            start = startdate, end = enddate, maxResults=500)
                            # flightDirection=ascdesc,
    if outAspd:
        # df=pd.DataFrame([p.properties for p in results])
        df=pd.DataFrame([{**p.properties, "geometry": shape(p.geometry)} for p in results])
        if df.empty:
            print('ASF returned empty output')
            return df
        cols = ['flightDirection','pathNumber','frameNumber','startTime','sceneName','url', 'geometry']
        df = df[cols]
        df=gpd.GeoDataFrame(df, geometry="geometry", crs="EPSG:4326")
        # download size is in pd.DataFrame.from_dict(r.bytes)
        # the url CAN be used directly with wget_alaska approach...
        return df
    else:
        return results

'''
# first get some results e.g. as
lat=33.74; lon=-118.37
wkt = f"POINT({lon} {lat})"
results=asf.geo_search(dataset='NISAR', processingLevel='GSLC', intersectsWith=wkt, maxResults=50)

# then you can download e.g. all of them, as:
images = []
downloadit=True
asf_session = get_asf_session()
downpath = 'NISAR'
reslen = len(results)
i = 0
for gg in results:
    i = i+1
    imname = gg.properties['sceneName']
    images.append(imname)
    if downloadit:
        print('Downloading '+imname+' ({0}/{1})'.format(str(i), str(reslen)))
        start = time.time()
        gg.download(downpath, session=asf_session)
        end = time.time()
        print(f"downloaded in {(end - start)/60:.2f} minutes")
'''


def list_sizes(path = 'NISAR_L2_PR_GSLC_009_034_A_018_4005_DHDH_A_20251230T130752_20251230T130827_X05009_N_F_J_001.h5'):
    ''' useful function to list contents of the H5 file - careful, seems the listed sizes are not exactly correct..
    '''
    with h5py.File(path, "r") as f:
        def visit(name, obj):
            if isinstance(obj, h5py.Dataset):
                size_bytes = obj.size * obj.dtype.itemsize
                size_mb = size_bytes / (1024**2)
                print(f"{name}: {size_mb:.3f} MB")
        f.visititems(visit)


def load_gslc(path, freq_code = 'A', polarization = 'HH', chunks="auto"):
    """
    Lazily load OPERA/GSLC FrequencyA HH grid as xarray.Dataset with:
    - complex data
    - x/y coordinates
    - EPSG CRS (if present)
    """
    f = h5py.File(path, "r")
    basestr = '/science/LSAR/GSLC/grids/frequency'+freq_code
    # +'/'+polarization
    # could do also from RSLC, but then the geocoding etc......
    # base='/science/LSAR/RSLC/swaths/frequencyA/HH'
    # --- main complex grid (lazy) ---
    dset = f[basestr+'/'+polarization]
    data = da.from_array(dset, chunks=chunks)
    # --- coordinates (lazy) ---
    x = da.from_array(f[basestr+"/xCoordinates"], chunks=chunks)
    y = da.from_array(f[basestr+"/yCoordinates"], chunks=chunks)
    # --- CRS ---
    proj_group = f[basestr+"/projection"]
    epsg = proj_group.attrs.get("epsg_code", None)
    crs = CRS.from_epsg(int(epsg)).to_string() if epsg is not None else None
    # --- Build xarray Dataset ---
    ds = xr.Dataset(
        data_vars={
            str(polarization): (("y", "x"), data)
        },
        coords={
            "x": ("x", x),
            "y": ("y", y),
        },
        attrs={
            "epsg": epsg,
            "crs": crs,
            "source_file": path,
            "freq": freq_code,
            "polarization": polarization
        }
    )
    return ds


'''
in1='NISAR_L2_PR_GSLC_009_034_A_018_4005_DHDH_A_20251230T130752_20251230T130827_X05009_N_F_J_001.h5'
ds1 = load_gslc(in1, freq_code = 'A')
in2='NISAR_L2_PR_GSLC_007_034_A_018_4005_DHDH_A_20251206T130751_20251206T130826_X05009_N_F_J_001.h5'
ds2 = load_gslc(in2, freq_code = 'A')

# /science/LSAR/GSLC/metadata/sourceData/swaths/frequencyA/nearRangeIncidenceAngle
# /science/LSAR/GSLC/metadata/sourceData/processingInformation/parameters/frequencyA/slantRange
# /science/LSAR/GSLC/metadata/sourceData/processingInformation/parameters/slantRange
# /science/LSAR/GSLC/metadata/sourceData/processingInformation/parameters/referenceTerrainHeight
# /science/LSAR/GSLC/metadata/radarGrid/incidenceAngle
science/LSAR/GSLC/metadata/radarGrid/losUnitVectorX
science/LSAR/GSLC/metadata/radarGrid/losUnitVectorY
science/LSAR/GSLC/metadata/radarGrid/elevationAngle
science/LSAR/GSLC/metadata/radarGrid/alongTrackUnitVectorX
science/LSAR/GSLC/metadata/radarGrid/alongTrackUnitVectorY
science/LSAR/GSLC/metadata/radarGrid/projection
science/LSAR/GSLC/metadata/radarGrid/xCoordinates
science/LSAR/GSLC/metadata/radarGrid/yCoordinates
# science/LSAR/GSLC/metadata/radarGrid/zeroDopplerAzimuthTime
'''


def get_ENU(path, chunks='auto'):
    ''' extracts the ENU unit vectors from the GSLC H5 file
    'NISAR_L2_PR_GSLC_009_034_A_018_4005_DHDH_A_20251230T130752_20251230T130827_X05009_N_F_J_001.h5'
    NEEDS SOME CHECKS!!! is it same direction as expected by our LiCS definition of the look angle data??
    '''
    chunks = 'auto'
    f = h5py.File(path, "r")
    basestr = 'science/LSAR/GSLC/metadata/radarGrid' #/xCoordinates'
    inc = f[basestr+'/incidenceAngle']
    inc = da.from_array(inc, chunks=chunks)
    unit_x = f[basestr+'/losUnitVectorX']
    unit_x = da.from_array(unit_x, chunks=chunks)
    unit_y = f[basestr+'/losUnitVectorY']
    unit_y = da.from_array(unit_y, chunks=chunks)
    # --- the shape is 21,y,x... and values similar - using mean then... ---
    inc=inc.mean(axis=[0])
    unit_x=unit_x.mean(axis=[0])
    unit_y=unit_y.mean(axis=[0])
    # --- coordinates (lazy) ---
    x = da.from_array(f[basestr+"/xCoordinates"], chunks=chunks)
    y = da.from_array(f[basestr+"/yCoordinates"], chunks=chunks)
    # --- CRS ---
    proj_group = f[basestr+"/projection"]
    epsg = proj_group.attrs.get("epsg_code", None)
    crs = CRS.from_epsg(int(epsg)).to_string() if epsg is not None else None
    # --- Build xarray Dataset ---
    ds = xr.Dataset(
        data_vars={
            # "inc": (("y", "x"), inc),   # or just unit_z to be cos(inc)
            "unit_x": (("y", "x"), unit_x),
            "unit_y": (("y", "x"), unit_y),
            "unit_z": (("y", "x"), np.cos(np.deg2rad(inc)))
        },
        coords={
            "x": ("x", x),
            "y": ("y", y),
        },
        attrs={
            "epsg": epsg,
            "crs": crs,
            "source_file": path,
        }
    )
    return ds



def generate_ifg(in1='NISAR_L2_PR_GSLC_009_034_A_018_4005_DHDH_A_20251230T130752_20251230T130827_X05009_N_F_J_001.h5',
                in2='NISAR_L2_PR_GSLC_007_034_A_018_4005_DHDH_A_20251206T130751_20251206T130826_X05009_N_F_J_001.h5',
                freq_code = 'A', polarization = 'HH', clipping_bbox = None,
                target_resolution_m = 110, outncfile = 'ifg.A.nc', create_wgs84_previews = True):
    ''' Main code to create interferogram and coherence from two NISAR GSLC data

    clipping_bbox should be shapely.Polygon in WGS-84
    '''
    if (not os.path.exists(in1)) or (not os.path.exists(in2)):
        print('ERROR - some of the files do not exist')
        return False
    try:
        ds1 = load_gslc(in1, freq_code=freq_code, polarization = polarization)
    except:
        print('ERROR loading epoch1. maybe wrong polarization? trying VV')
        polarization = 'VV'
        ds1 = load_gslc(in1, freq_code=freq_code, polarization=polarization)
    ds2 = load_gslc(in2, freq_code=freq_code, polarization = polarization)
    if type(clipping_bbox) != type(None):
        utmcode = ds1.attrs.get("crs")
        clipping_bbox_utm = wgs2utm(clipping_bbox, utmcode)
        x1, y1, x2, y2 = clipping_bbox_utm.bounds
        ds1 = ds1.sel(x=slice(x1,x2),
                      y=slice(y2,y1))
        ds2 = ds2.sel(x=slice(x1, x2),
                      y=slice(y2, y1))
    if not ds1[polarization].shape == ds2[polarization].shape:
        print('ERROR - the dimensions do not fit - are these files from the same track and path??')
        return False
    # Lazy interferogram
    ifg = ds1.HH * ds2.HH.conj()
    #
    ifg_da = xr.DataArray(
        ifg,
        coords=ds1.coords,
        dims=ds1[polarization].dims,
        attrs={"description": "Interferogram (S1 * conj(S2))"}
    )
    #
    ifg_da.attrs.update({
        "source1": ds1.attrs.get("source_file", os.path.basename(in1)),
        "source2": ds2.attrs.get("source_file", os.path.basename(in2)),
        "crs": ds1.attrs.get("crs"),
        "epsg": ds1.attrs.get("epsg"),
        "freq": ds1.attrs.get("freq"),
        "polarization": ds1.attrs.get("polarization"),
    })
    #
    # now multilook
    # target_res = 110 # [m] first instance, to get similar outputs to LiCSAR ifgs..
    resx = float(ifg_da.x[1]-ifg_da.x[0])
    resy = abs(float(ifg_da.y[1]-ifg_da.y[0]))
    mlx=int(target_resolution_m/resx) # setting the floor here - may still be too large then?
    mly=int(target_resolution_m/resy)
    print(f'Using {mlx}x{mly} multilooking, i.e. coherence calculated based on {mlx*mly} pixels')
    ifg_ml_da = ifg_da.coarsen({'x': mlx, 'y': mly}, boundary='trim')
    
    # then, create ~coherence estimate (see lics_unwrap.py)
    amp = np.abs(ifg_da).coarsen({'x': mlx, 'y': mly}, boundary='trim')
    coh = np.abs(ifg_ml_da.sum()) / amp.sum()  # quality of this coherence depends on multilook pixels!

    # multilooked phase
    ifg_ml = ifg_ml_da.sum()/ifg_ml_da.count()
    phase = xr.apply_ufunc(
        np.angle,
        ifg_ml,
        dask="parallelized",
        output_dtypes=[np.float32],
    )

    # wrap it up to dataset
    ds_out = xr.Dataset(
        data_vars={
            "phase": phase, #(ifg_ml.dims, phase.data),
            "coherence": coh, # (ifg_da.dims, magnitude.data),
        },
        coords=ifg_ml.coords,
        attrs={
            "description": "Interferogram components from S1 * conj(S2)",
            "source1": str(ifg_ml.attrs.get("source1", "-")),
            "source2": str(ifg_ml.attrs.get("source2", "-")),
            "crs": ifg_ml.attrs.get("crs", '-'),
            "epsg": ifg_ml.attrs.get("epsg", '-'),
            "freq": ifg_ml.attrs.get("freq", '-'),
            "polarization": ifg_ml.attrs.get("polarization", '-'),
        },
    )

    chunksizes = tuple(ch[0] for ch in phase.data.chunks)  # e.g. (1024, 1024)

    encoding = {
        "phase": {
            "zlib": True, "complevel": 4,  # netCDF4/deflate compression
            "dtype": "float32",
            "chunksizes": chunksizes,       # ensure chunk-wise writing
            "_FillValue": np.float32(np.nan),
        },
        "coherence": {
            "zlib": True, "complevel": 4,
            "dtype": "float32",
            "chunksizes": chunksizes,
            "_FillValue": np.float32(np.nan),
        },
    }
    print('Computing and storing to '+outncfile)
    # This will execute lazily by chunks, not loading everything into memory
    ds_out.to_netcdf(outncfile, engine="netcdf4", encoding=encoding)
    if create_wgs84_previews:
        # this means first convert to WGS84 and then .. create some nice pngs
        # just in case of rio issues, reloading:
        xrds = xr.open_dataset(outncfile)
        convert_to_wgs84(xrds, outbasename = outncfile[:-3])
        # print('TODO - but see the code here..')
        print('all done')

'''
# Lazily compute phase and magnitude (no .values!)
phase = xr.apply_ufunc(
    np.angle,
    ifg_da,
    dask="parallelized",
    output_dtypes=[np.float32],
)

magnitude = xr.apply_ufunc(
    np.abs,
    ifg_da,
    dask="parallelized",
    output_dtypes=[np.float32],
)

ds_out = xr.Dataset(
    data_vars={
        "phase": (ifg_da.dims, phase.data),
        "magnitude": (ifg_da.dims, magnitude.data),
    },
    coords=ifg_da.coords,
    attrs={
        "description": "Interferogram components from S1 * conj(S2)",
        #"source1": str(ifg_da.attrs.get("source1", "")),
        #"source2": str(ifg_da.attrs.get("source2", "")),
        "crs": ifg_da.attrs.get("crs", None),
        "epsg": ifg_da.attrs.get("epsg", None),
    },
)

# Choose compression + chunk sizes for NetCDF
# Use the same chunking as in the dask graph to avoid rechunking during write
# xarray expects chunk sizes per dimension as a tuple
chunksizes = tuple(ch[0] for ch in phase.data.chunks)  # e.g. (1024, 1024)

encoding = {
    "phase": {
        "zlib": True, "complevel": 4,  # netCDF4/deflate compression
        "dtype": "float32",
        "chunksizes": chunksizes,       # ensure chunk-wise writing
        "_FillValue": np.float32(np.nan),
    },
    "magnitude": {
        "zlib": True, "complevel": 4,
        "dtype": "float32",
        "chunksizes": chunksizes,
        "_FillValue": np.float32(np.nan),
    },
}

# This will execute lazily by chunks, not loading everything into memory
ds_out.to_netcdf("ifg_components.nc", engine="netcdf4", encoding=encoding)

'''

# but the orig data is still large - downsample to 50x50 m:
'''
def multilook_ifg_phase_mag_to_netcdf(
    nc_in: str,
    nc_out: str,
    phase_var: str = "phase",
    magnitude_var: str = "magnitude",
    looks_y: int = 10,
    looks_x: int = 10,
    chunks: dict | None = None,
    write_magnitude: bool = True,
    netcdf_engine: str = "netcdf4",
    coord_reduce: str = "mean"):
    """
    Lazily multilook an interferogram stored as phase/magnitude in NetCDF and
    write the multilooked phase (and optionally magnitude) to a new NetCDF.

    Strategy:
        1) Reconstruct complex: C = M * exp(i * phase)
        2) Coarsen (mean) with factors (looks_y, looks_x), boundary='trim'
        3) Phase_out = angle(mean(C)), Magnitude_out = abs(mean(C))
        4) Coarsen 1D coords x, y to match new sizes (mean or window center)
        5) Write to NetCDF with chunked encoding (out-of-core)

    Notes:
        - Assumes 2D variables with dims ("y", "x") or a permutation thereof.
        - Coordinates 'x' and 'y' are assumed to be 1D monotonic vectors.
    """
    # 1) Open lazily
    ds = xr.open_dataset(nc_in, engine=netcdf_engine, chunks=chunks)

    # Identify dims of the phase variable (e.g., ('y','x') or ('x','y'))
    pdims = ds[phase_var].dims
    if len(pdims) != 2:
        raise ValueError(f"{phase_var} must be 2D, got dims {pdims}")

    # Make sure order is (y, x)
    if pdims == ("y", "x"):
        ydim, xdim = "y", "x"
        phase = ds[phase_var]
        mag   = ds[magnitude_var]
    elif pdims == ("x", "y"):
        ydim, xdim = "y", "x"
        phase = ds[phase_var].transpose("y", "x")
        mag   = ds[magnitude_var].transpose("y", "x")
    else:
        # Fall back: sort dims by name if not exact
        dims_sorted = tuple(sorted(pdims))
        raise ValueError(f"Unexpected dims {pdims}; expected ('y','x') or ('x','y').")

    phase = phase.astype("float32")
    mag   = mag.astype("float32")

    # 2) Reconstruct complex interferogram lazily: C = M * exp(i * phase)
    C = mag * xr.apply_ufunc(
        lambda p: np.exp(1j * p),
        phase,
        dask="parallelized",
        output_dtypes=[np.complex64],
    ).astype("complex64")

    # 3) Multilook via coarsen (mean), trimming edges that don't complete a window
    C_ml = C.coarsen({ydim: looks_y, xdim: looks_x}, boundary="trim").mean()

    # Convert back to phase and magnitude lazily
    phase_ml = xr.apply_ufunc(np.angle, C_ml, dask="parallelized", output_dtypes=[np.float32]).astype("float32")
    if write_magnitude:
        magnitude_ml = xr.apply_ufunc(np.abs, C_ml, dask="parallelized", output_dtypes=[np.float32]).astype("float32")

    # 4) Build coarsened coordinates that MATCH the multilooked sizes
    coords_out = {}

    # Coarsen y coordinate if present & 1D
    if ydim in ds.coords and ds[ydim].ndim == 1:
        if coord_reduce == "mean":
            y_ml = ds[ydim].coarsen({ydim: looks_y}, boundary="trim").mean().astype("float32")
        elif coord_reduce == "center":
            n_y = ds.sizes[ydim]
            n_y_trim = (n_y // looks_y) * looks_y
            y_trim = ds[ydim].isel({ydim: slice(0, n_y_trim)})
            y_ml = y_trim.coarsen({ydim: looks_y}).construct(f"{ydim}_win").isel({f"{ydim}_win": looks_y // 2})
            y_ml = y_ml.astype("float32")
        else:
            raise ValueError("coord_reduce must be 'mean' or 'center'")
        coords_out[ydim] = (ydim, y_ml.data)

    # Coarsen x coordinate if present & 1D
    if xdim in ds.coords and ds[xdim].ndim == 1:
        if coord_reduce == "mean":
            x_ml = ds[xdim].coarsen({xdim: looks_x}, boundary="trim").mean().astype("float32")
        elif coord_reduce == "center":
            n_x = ds.sizes[xdim]
            n_x_trim = (n_x // looks_x) * looks_x
            x_trim = ds[xdim].isel({xdim: slice(0, n_x_trim)})
            x_ml = x_trim.coarsen({xdim: looks_x}).construct(f"{xdim}_win").isel({f"{xdim}_win": looks_x // 2})
            x_ml = x_ml.astype("float32")
        else:
            raise ValueError("coord_reduce must be 'mean' or 'center'")
        coords_out[xdim] = (xdim, x_ml.data)

    # 5) Output dataset with matching dims & coords
    data_vars = {"phase_ml": (phase_ml.dims, phase_ml.data)}
    if write_magnitude:
        data_vars["magnitude_ml"] = (phase_ml.dims, magnitude_ml.data)

    # Safe attrs (avoid None)
    attrs = {"description": f"Multilooked {looks_y}x{looks_x} from complex mean"}
    for key in ("crs", "epsg", "source1", "source2"):
        val = ds.attrs.get(key)
        if val is not None:
            attrs[key] = str(val)

    ds_out = xr.Dataset(data_vars=data_vars, coords=coords_out, attrs=attrs)

    # 6) Chunked encoding (use first chunk per dimension)
    out_chunks = tuple(ch[0] for ch in phase_ml.data.chunks)  # e.g., (1024, 1024)
    encoding = {
        "phase_ml": {
            "zlib": True, "complevel": 4,
            "dtype": "float32",
            "chunksizes": out_chunks,
            "_FillValue": np.float32(np.nan),
        }
    }
    if write_magnitude:
        encoding["magnitude_ml"] = {
            "zlib": True, "complevel": 4,
            "dtype": "float32",
            "chunksizes": out_chunks,
            "_FillValue": np.float32(np.nan),
        }

    # 7) Write lazily, chunk-by-chunk (no full-RAM load)
    ds_out.to_netcdf(nc_out, engine=netcdf_engine, encoding=encoding)
    print('stored to '+nc_out)


# Example: 10×10 multilook, chunks are multiples of 10 for efficiency
multilook_ifg_phase_mag_to_netcdf(
    nc_in="ifg_components.nc",
    nc_out="ifg_components_ml_10x10.nc",
    looks_y=10,
    looks_x=10,
    chunks={"y": 4000, "x": 4000},  # multiples of looks are efficient
    write_magnitude=True,
    coord_reduce="mean",            # or "center"
)
'''

# now convert to wgs-84 geotiff:
def convert_to_wgs84(xrds, outbasename = 'ifg', create_previews = True):
    ''' This will convert the xr.Dataset containing phase and coherence to WGS-84 geotiffs
    you can either parse the ds_out, or use e.g. xrds=xr.open_dataset('ifg.nc') 
    '''
    if 'coherence' not in xrds:
        print('ERROR: coherence is not part of the input data (is this xr.Dataset?)')
        return False
    p=xrds.phase
    p=p.rio.write_crs(xrds.crs)
    p.rio.to_raster(outbasename+'_pha.tif')
    m=xrds.coherence #magnitude_ml
    m=m.rio.write_crs(xrds.crs)
    m.rio.to_raster(outbasename+'_coh.tif') # no need for compression as i will translate to wgs later
    cmd = "gdalwarp -t_srs EPSG:4326 -r near -co COMPRESS=DEFLATE -co PREDICTOR=2 "+outbasename+"_pha.tif "+outbasename+"_pha.wgs84.tif"
    os.system(cmd)
    cmd = "gdalwarp -t_srs EPSG:4326 -r average -co COMPRESS=DEFLATE -co PREDICTOR=2 "+outbasename+"_coh.tif "+outbasename+"_coh.wgs84.tif"
    os.system(cmd)
    if create_previews:
        # now you can create e.g. some nice previews
        cmd = "create_preview_pygmt.py --grid "+outbasename+"_pha.wgs84.tif --title NISAR --label phase --photobg --lims -3.141593 3.141593"
        os.system(cmd)
        cmd = "create_preview_pygmt.py --grid "+outbasename+"_coh.wgs84.tif --title NISAR --cmap gray --label coherence --photobg --lims 0 1"
        os.system(cmd)


# this below is TODO for the unit vectors:
'''
z=a.unit_z
z=z.rio.write_crs(a.crs)
z.rio.to_raster('unitz.tif')
cmd = "gdalwarp -t_srs EPSG:4326 -r near -co COMPRESS=DEFLATE -co PREDICTOR=2 unitz.tif unitz.wgs84.tif"
os.system(cmd)
'''

'''
# then to convert to WGS-84:
gdal_translate \
  -of GTiff \
  -co COMPRESS=LZW \
  -co TILED=YES \
  NETCDF:ifg_components.nc:phase \
  phase_utm.tif

gdalwarp \
  -t_srs EPSG:4326 \
  -r nearest \
  -co COMPRESS=DEFLATE -co PREDICTOR=2 \
  -co TILED=YES \
  phase_utm.tif \
  phase_wgs84.tif
'''


def download(url, slcdir = '/gws/ssde/j25a/nceo_geohazards/vol2/LiCS/temp/SLC', provider='alaska'):
    '''wrapper to wget commands. the provider must be one of ['cdse', 'alaska']
    '''
    # slcdir = os.environ['LiCSAR_SLC']
    wgetpath = os.environ['LiCSARpath']+'/bin/scripts/wget_asf_nisar'
    cmd = 'cd {0}; {1} {2}'.format(slcdir, wgetpath, url)
    rc = os.system(cmd)
    filepath = os.path.join(slcdir,filename)
    return filepath


def get_asf_session():
    credentials = os.path.join(os.environ['HOME'], '.asf_credentials')
    if not os.path.exists(credentials):
        credentials = os.path.join(os.environ["LiCSAR_configpath"],'asf_credentials')
    f = open(credentials,'r')
    asf_user = re.sub(r'\W+','',f.readline().split('=')[1])
    asf_pass = re.sub(r'\W+','',f.readline().split('=')[1])
    f.close()
    return asf.ASFSession().auth_with_creds(asf_user, asf_pass)

