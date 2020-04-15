#!/usr/bin/env python3
"""
========
Overview
========
This script converts GACOS products (.ztd and .ztd.rsc) to GeoTIFF files of slantrange and zenith total delay data.

v1.1 20200228 Yu Morishita, Uni of Leeds and GSI

=====
Usage
=====
LiCSAR_ztd2geotiff.py -z ztdfile -u geo.U.tif 

 -z  Path to input ztd file (float, little endian, with .ztd.rsc file)
 -u  Path to GeoTIFF file of vertical component of LOS unit vector

"""

#%% Change log
'''
v1.1 20200228 Yu Morishita, Uni of Leeds and GSI
 - Output png
v1.0 20200218 Yu Morishita, Uni of Leeds and GSI
 - Original implementationf
'''

#%% Import
import getopt
import os
import sys
import time
import subprocess as subp
import numpy as np
import gdal, osr

os.environ['QT_QPA_PLATFORM']='offscreen'
import warnings
import matplotlib
with warnings.catch_warnings(): ## To silence user warning
    warnings.simplefilter('ignore', UserWarning)
    matplotlib.use('Agg')
from matplotlib import pyplot as plt

class Usage(Exception):
    """Usage context manager"""
    def __init__(self, msg):
        self.msg = msg

#%%
def get_param_par(mlipar, field):
    value = subp.check_output(['grep', field,mlipar]).decode().split()[1].strip()
    return value

#%%
def make_geotiff(data, length, width, latn, lonw, dlat, dlon, outfile, compress_option):
    nodata = np.nan  ## or 0?
    
    ## Grid registration to pixel registration by shifing half pixel
    latn_p = latn - dlat/2
    lonw_p = lonw - dlon/2

    driver = gdal.GetDriverByName('GTiff')
    outRaster = driver.Create(outfile, width, length, 1, gdal.GDT_Float32, options=compress_option)
    outRaster.SetGeoTransform((lonw_p, dlon, 0, latn_p, 0, dlat))
    outband = outRaster.GetRasterBand(1)
    outband.WriteArray(data)
    outband.SetNoDataValue(nodata)
    outRaster.SetMetadataItem('AREA_OR_POINT', 'Point')
    outRasterSRS = osr.SpatialReference()
    outRasterSRS.ImportFromEPSG(4326)
    outRaster.SetProjection(outRasterSRS.ExportToWkt())
    outband.FlushCache()

    return

#%%
def make_png(data, pngfile, cmap, title, vmin=None, vmax=None, cbar=True):
    """
    Make png image.
    """
    
    length, width = data.shape
    figsizex = 8
    xmergin = 2 if cbar else 0
    figsizey = int((figsizex-xmergin)*(length/width))+1
    
    ### Plot
    fig, ax = plt.subplots(1, 1, figsize=(figsizex, figsizey))
    plt.tight_layout()
    
    im = ax.imshow(data, vmin=vmin, vmax=vmax, cmap=cmap)
    ax.set_xticklabels([])
    ax.set_yticklabels([])
    ax.set_title(title)
    if cbar: fig.colorbar(im)

    plt.savefig(pngfile)
    plt.close()
    
    return

#%% Main
def main(argv=None):
        
    #%% Check argv
    if argv == None:
        argv = sys.argv
        
    start = time.time()
    ver=1.1; date=20200228; author="Y. Morishita"
    print("\n{} ver{} {} {}".format(os.path.basename(argv[0]), ver, date, author), flush=True)
    print("{} {}".format(os.path.basename(argv[0]), ' '.join(argv[1:])), flush=True)

    #%% Set default
    inztdfile = []
    LOSufile = []

    compress_option = ['COMPRESS=DEFLATE', 'PREDICTOR=3']
    ## ['COMPRESS=LZW', 'PREDICTOR=3'], ['COMPRESS=PACKBITS']
    resampleAlg = 'cubicspline'# None # 'cubic' 
    radar_frequency = 5405000000 #Hz, Sentinel-1
    speed_of_light = 299792458 #m/s
    wavelength = speed_of_light/radar_frequency #meter
    m2r_coef = 4*np.pi/wavelength
    cmap = 'inferno'
    

    #%% Read options
    try:
        try:
            opts, args = getopt.getopt(argv[1:], "hz:u:", ["help"])
        except getopt.error as msg:
            raise Usage(msg)
        for o, a in opts:
            if o == '-h' or o == '--help':
                print(__doc__)
                return 0
            elif o == '-z':
                inztdfile = a
            elif o == '-u':
                LOSufile = a

        rscfile = inztdfile+'.rsc'
        if not inztdfile:
            raise Usage('No input ztd file given, -z is not optional!')
        elif not os.path.exists(inztdfile):
            raise Usage('No {} exists!'.format(inztdfile))
        elif not os.path.exists(rscfile):
            raise Usage('No {} exists!'.format(rscfile))
        elif not LOSufile:
            raise Usage('No LOSu file given, -u is not optional!')
        elif not os.path.exists(LOSufile):
            raise Usage('No {} exists!'.format(LOSufile))

    except Usage as err:
        print("\nERROR:", file=sys.stderr, end='')
        print("  "+str(err.msg), file=sys.stderr)
        print("\nFor help, use -h or --help.\n", file=sys.stderr)
        return 2


    #%% Read LOSu file
    geotiff = gdal.Open(LOSufile)
    LOSu = geotiff.ReadAsArray()
    width_geo = geotiff.RasterXSize
    length_geo = geotiff.RasterYSize
    lonw_p, dlon_geo, _, latn_p, _, dlat_geo = geotiff.GetGeoTransform()
    ## lat lon are in pixel registration. dlat is negative
    latn_geo = latn_p + dlat_geo/2
    lonw_geo = lonw_p + dlon_geo/2
    lats_geo = latn_geo+dlat_geo*(length_geo-1)
    lone_geo = lonw_geo+dlon_geo*(width_geo-1)

    LOSu[LOSu==0] = np.nan


    #%% Read ztd. Grid registration
    width_ztd = int(get_param_par(rscfile, 'WIDTH'))
    length_ztd = int(get_param_par(rscfile, 'FILE_LENGTH'))
    dlat_ztd = float(get_param_par(rscfile, 'Y_STEP')) #minus
    dlon_ztd = float(get_param_par(rscfile, 'X_STEP'))
    latn_ztd = float(get_param_par(rscfile, 'Y_FIRST'))
    lonw_ztd = float(get_param_par(rscfile, 'X_FIRST'))


    #%% Make hdr file of ztd
    strings = ["NROWS          {}".format(length_ztd),
               "NCOLS          {}".format(width_ztd),
               "NBITS          32",
               "PIXELTYPE      FLOAT",
               "BYTEORDER      I",
               "LAYOUT         BIL",
               "ULXMAP         {}".format(lonw_ztd),
               "ULYMAP         {}".format(latn_ztd),
               "XDIM           {}".format(dlon_ztd),
               "YDIM           {}".format(np.abs(dlat_ztd))]
    hdrfile = inztdfile+'.hdr'
    with open(hdrfile, "w") as f:
        f.write("\n".join(strings))


    #%% Process ztd files 
    bilfile = inztdfile+'.bil'
    if os.path.exists(bilfile): os.remove(bilfile)
    os.symlink(inztdfile, bilfile)

    ### Cut and resapmle ztd to geo
    print('\nCut and resapmle ztd...', flush=True)
    ztd_geo = gdal.Warp("", bilfile, format='MEM', outputBounds=(lonw_geo, lats_geo, lone_geo, latn_geo), width=width_geo, height=length_geo, resampleAlg=resampleAlg, srcNodata=0).ReadAsArray()
    ztd_geo = ztd_geo*m2r_coef ## meter -> rad
    #ztd_geo[ztd_geo==0] = np.nan
    
    print('\nCompute sltd from ztd and LOSu...', flush=True)
    sltd_geo = ztd_geo/LOSu ## LOSu=cos(inc)

    os.remove(hdrfile)
    os.remove(bilfile)


    #%% Output geotiff
    print('\nOutput GeoTIFF files...', flush=True)
    outztdfile = inztdfile+'.geo.tif'
    outsltdfile = outztdfile.replace('ztd', 'sltd')

    make_geotiff(ztd_geo, length_geo, width_geo, latn_geo, lonw_geo, dlat_geo, dlon_geo, outztdfile, compress_option)
    make_geotiff(sltd_geo, length_geo, width_geo, latn_geo, lonw_geo, dlat_geo, dlon_geo, outsltdfile, compress_option)

    
    #%% Create png
    ztdpng = inztdfile+'.geo.png'
    sltdpng = ztdpng.replace('ztd', 'sltd')

    make_png(ztd_geo, ztdpng, cmap, 'ztd (rad)', vmin=None, vmax=None, cbar=True)
    make_png(sltd_geo, sltdpng, cmap, 'sltd (rad)', vmin=None, vmax=None, cbar=True)


    #%% Finish
    elapsed_time = time.time()-start
    hour = int(elapsed_time/3600)
    minite = int(np.mod((elapsed_time/60),60))
    sec = int(np.mod(elapsed_time,60))
    print("\nElapsed time: {0:02}h {1:02}m {2:02}s".format(hour,minite,sec))

    print('\n{} Successfully finished!!\n'.format(os.path.basename(argv[0])))
    print('Output: {}, {}\n'.format(outztdfile, outsltdfile), flush=True)
    print('        {}, {}\n'.format(ztdpng, sltdpng), flush=True)


#%% main
if __name__ == "__main__":
    sys.exit(main())
