#!/usr/bin/env python3
'Muhammet Nergizci, COMET University of Leeds, 2025'
import argparse
import os
import time
import sys
import rioxarray
import numpy as np
import xarray as xr

def parse_args():
    parser = argparse.ArgumentParser(
        description="This script prepares the ENU vector in the azimuth direction."
    )
    parser.add_argument("frame", help="The frame name from the LiCS portal, e.g., 094D_05288_130913")
    return parser.parse_args()


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
    xrpha = xrpha.drop_vars('band')
    
    if cliparea_geo:
        minclipx, maxclipx, minclipy, maxclipy = cliparea_geo2coords(cliparea_geo)
        xrpha = xrpha.sel(x=slice(minclipx, maxclipx), y=slice(maxclipy, minclipy))
    if tolonlat:
        xrpha = xrpha.rename({'x': 'lon','y': 'lat'})
    return xrpha

def export_xr2tif(xrda, tif, lonlat = True, debug = True, dogdal = True, refto = None, set_to_pixel_registration = False):
    """Exports xarray dataarray to a geotiff
    
     Args:
        xrda (xarray.Dataarray): dataarray to export
        tif (string): path to output tif file
        lonlat (boolean): are the dimensions named as lon, lat?
        debug (boolean): just load it as float32
        dogdal (boolean): after exporting, perform gdalwarp (fix for potential issues in output geotiff) and gdal_translate to better compress
        refto (str): path to the (usually hgt) file to apply gdalwarp2match.py to. If None, it will not apply
        set_to_pixel_registration (boolean):  will rewrite header to assume Pixel Registration - that's by default in LiCSAR data (but not in rasterio...)
    """
    import rioxarray
    #coordsys = xrda.crs.split('=')[1]
    coordsys = "epsg:4326"
    if debug:
        xrda = xrda.astype(np.float32)
        # reset original spatial_ref
        if 'spatial_ref' in xrda:
            xrda = xrda.drop_vars('spatial_ref')
        # remove attributes
        xrda.attrs = {}
    if lonlat:
        if xrda.lat[1] > xrda.lat[0]:
            xrda = xrda.sortby('lat',ascending=False)
        xrda = xrda.rio.set_spatial_dims(x_dim="lon", y_dim="lat", inplace=True)
    else:
        if xrda.y[1] > xrda.y[0]:
            xrda = xrda.sortby('y',ascending=False)
        xrda = xrda.rio.set_spatial_dims(x_dim="x", y_dim="y", inplace=True)
    xrda = xrda.rio.write_crs(coordsys, inplace=True)
    if dogdal:
        xrda.rio.to_raster(tif+'tmp.tif')
        if refto:
            cmd = 'gdalwarp2match.py {0} {1} {2}; mv {2} {0}'.format(tif+'tmp.tif', refto, tif)
            runcmd(cmd, printcmd=False)
        else:
            cmd = 'gdalwarp -t_srs EPSG:4326 {0} {1}'.format(tif + 'tmp.tif', tif)  # will fix some issues
            runcmd(cmd, printcmd=False)
        if set_to_pixel_registration:
            cmd = 'gdal_edit.py -mo AREA_OR_POINT=Point '+tif
            runcmd(cmd, printcmd=False)
        cmd = 'mv {0} {1}; gdal_translate -of GTiff -co COMPRESS=DEFLATE -co PREDICTOR=3 {1} {0}'.format(tif, tif+'tmp.tif') # will compress better
        runcmd(cmd, printcmd = False)
        if os.path.exists(tif+'tmp.tif'):
            os.remove(tif+'tmp.tif')
        if os.path.exists(tif+'tmp2.tif'):
            os.remove(tif + 'tmp2.tif')
    else:
        xrda.rio.to_raster(tif, compress='deflate')
        if set_to_pixel_registration:
            cmd = 'gdal_edit.py -mo AREA_OR_POINT=Point '+tif
            runcmd(cmd, printcmd=False)

def main():
    start_time = time.time()

    args = parse_args()
    frame_id = args.frame
    track = int(frame_id[0:3])
    orientation = frame_id[3]

    print(orientation)

    homedir = os.getcwd()
    LiCS_public = os.environ['LiCSAR_public']

    frame_dir = os.path.join(homedir, frame_id)
    geoc_dir = os.path.join(frame_dir, 'GEOC')

    print(f"Starting script for LiCSBAS frame: {frame_dir}")

    if not os.path.isdir(frame_dir):
        print(f"Frame directory '{frame_dir}' does not exist, creating..")
        os.makedirs(frame_dir, exist_ok=True)

    if not os.path.isdir(geoc_dir):
        print(f"ERROR: '{geoc_dir}' does not exist! The ENU will be taken from LiCS portal", file=sys.stderr)

        E_tif = os.path.join(LiCS_public, str(track), frame_id, 'metadata', f'{frame_id}.geo.E.tif')
        N_tif = os.path.join(LiCS_public, str(track), frame_id, 'metadata', f'{frame_id}.geo.N.tif')
        U_tif = os.path.join(LiCS_public, str(track), frame_id, 'metadata', f'{frame_id}.geo.U.tif')
    else:
        E_tif = os.path.join(geoc_dir, f'{frame_id}.geo.E.tif')
        N_tif = os.path.join(geoc_dir, f'{frame_id}.geo.N.tif')
        U_tif = os.path.join(geoc_dir, f'{frame_id}.geo.U.tif')

    print(E_tif)
    print(N_tif)
    print(U_tif)

    E = load_tif2xr(E_tif)
    N = load_tif2xr(N_tif)
    U = load_tif2xr(U_tif)

    E = xr.where(E == 0, np.nan, E)
    N = xr.where(N == 0, np.nan, N)
    U = xr.where(U == 0, np.nan, U)

    if orientation == "A":
        inc_rad = np.arccos(U)
        head_rad = np.arcsin(N / np.sin(inc_rad))
    elif orientation == "D":
        inc_rad = np.arccos(U)
        head_rad = np.arcsin(-N / np.sin(inc_rad)) - np.pi
    else:
        raise ValueError("The 4th character of frameID is neither A nor D, please check your frame name.")

    azi_N = np.cos(head_rad)
    azi_E = np.sin(head_rad)
    azi_U = xr.zeros_like(azi_E)

    if not os.path.isdir(geoc_dir):
        os.makedirs(geoc_dir, exist_ok=True)

    export_xr2tif(azi_E, os.path.join(geoc_dir, f"{frame_id}.geo.E.azi.tif"))
    export_xr2tif(azi_N, os.path.join(geoc_dir, f"{frame_id}.geo.N.azi.tif"))
    export_xr2tif(azi_U, os.path.join(geoc_dir, f"{frame_id}.geo.U.azi.tif"))

    meta_dir = os.path.join(LiCS_public, str(track), frame_id, 'metadata')
    if not os.path.exists(os.path.join(meta_dir, f"{frame_id}.geo.U.azi.tif")):
        export_xr2tif(azi_E, os.path.join(meta_dir, f"{frame_id}.geo.E.azi.tif"))
        export_xr2tif(azi_N, os.path.join(meta_dir, f"{frame_id}.geo.N.azi.tif"))
        export_xr2tif(azi_U, os.path.join(meta_dir, f"{frame_id}.geo.U.azi.tif"))

    print("Processing complete.")
    print(f"Elapsed time: {time.time() - start_time:.2f} seconds")
    return 0


if __name__ == "__main__":
    sys.exit(main())