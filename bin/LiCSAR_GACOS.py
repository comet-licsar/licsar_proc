#!/usr/bin/env python3
"""
========
Overview
========
This script applies a GACOS correction to an unw or diff data. GACOS data must be prepared beforehand by requesting on the GACOS web.

=========
Changelog
=========
v1.0 20191016 Yu Morishita, Uni of Leeds and GSI
 - Original implementation

=====
Usage
=====
LiCSAR_GACOS.py [-u unwfile | -d diffile] -m ztdfile_master -s ztdfile_slave -l [geo.U.tif | incidence_angle] [-p EQA.dem_par] [--bigendian]

 -u  Input unw file to be corrected (must be geocoded).
     Format can be GeoTIFF or float32 (phase only).
     If float32, -p option must be used to specify geographical coordinates.
 -d  Input diff file to be corrected (must be geocoded).
     Format can be GeoTIFF (phase only), complex64 (amp+phase),
     or float32 (phase only).
     If complex64 or float32, -p option must be used to specify
     geographical coordinates.
     Either -u or -d must be given.
 -m  Input ztd file (GACOS data) for master.
     Accompanying rsc file is also required.
 -s  Input ztd file (GACOS data) for slave.
 -l  GeoTIFF or float32 file of vertical component of LOS vector
     or incidence angle (scholar in degree).
 -p  EQA.dem_par file containing geographical coordinates.
     Use if input unw/diff file is in float32 or complex64.
 --bigendian  Use if unw/diff file is in float32/complex64 and bigendian.

Exapmles:
1) LiCSAR_GACOS.py -u 20170628_20170704.geo.unw.tif -m 20170628.ztd -s 20170704.ztd -l 014A_05138_131313.geo.U.tif
2) LiCSAR_GACOS.py -u 20170628_20170704.geo.unw -m 20170628.ztd -s 20170704.ztd -l 40 -p EQA.dem_par --bigendian
3) LiCSAR_GACOS.py -d 20170628_20170704.geo.diff_phase.tif -m 20170628.ztd -s 20170704.ztd -l U.geo
4) LiCSAR_GACOS.py -d 20170628_20170704.geo.diff -m 20170628.ztd -s 20170704.ztd -l 014A_05138_131313.geo.U.tif -p EQA.dem_par --bigendian
5) LiCSAR_GACOS.py -d 20170628_20170704.geo.diff_phase -m 20170628.ztd -s 20170704.ztd -l 014A_05138_131313.geo.U.tif -p EQA.dem_par --bigendian

"""

#%% Import
import getopt
import os
os.environ['QT_QPA_PLATFORM']='offscreen'
import sys
import time
import shutil
import numpy as np
import gdal, osr
import subprocess as subp

import matplotlib
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
def read_img(file, length, width, dtype=np.float32, endian='little'):
    if endian == 'little':
        data = np.fromfile(file, dtype=dtype).reshape((length, width))
    else:
        data = np.fromfile(file, dtype=dtype).byteswap().reshape((length, width))
    return data

#%%
def write_img(data, outfile, endian='little'):
    if endian == 'little':
        data.tofile(outfile)
    else:
        data.byteswap().tofile(outfile)
    return

#%%
def make_3im_png(data3, pngfile, cmap, title3, vmin=None, vmax=None, cbar=True):
    """
    Make png with 3 images for comparison.
    data3 and title3 must be list with 3 elements.
    cmap can be 'insar'. To wrap data, np.angle(np.exp(1j*x/cycle)*cycle)
    """
    ### Plot setting
    if cmap=='insar':
        cdict = cmap_insar()
        plt.register_cmap(name='insar',data=cdict)

    length, width = data3[0].shape
    figsizex = 12
    xmergin = 4 if cbar else 0
    figsizey = int((figsizex-xmergin)/3*length/width)+1

    fig = plt.figure(figsize = (figsizex, figsizey))

    for i in range(3):
        ax = fig.add_subplot(1, 3, i+1) #index start from 1
        im = ax.imshow(data3[i], vmin=vmin, vmax=vmax, cmap=cmap)
        ax.set_title(title3[i])
        ax.set_xticklabels([])
        ax.set_yticklabels([])
        if cbar: fig.colorbar(im, ax=ax)

    plt.tight_layout()
    plt.savefig(pngfile)
    plt.close()

    return
    
#%%
def cmap_insar():
    """
    How to use cmap_insar:
        cdict = cmap_insar()
        plt.register_cmap(name='insar', data=cdict)
        plt.imshow(array, cmap='insar', vmin=-np.pi, vmax=np.pi)
        
    Note:
        - Input array should be wrapped and in radian
        - To wrap unwrapped phase, np.angle(np.exp(1j*unw/cycle)*cycle)
    """

#    These are for 0-2pi, not -pi to pi
    red = [255,255,255,255,240,203,165,127,90,55,92,130,167,205,243,255,255,255]
    green = [118,156,193,231,255,255,255,255,255,255,217,179,142,104,66,80,118,156]
    blue = [191,153,116,78,69,106,144,182,219,255,255,255,255,255,255,229,191,153]
    phase = [k/32 for k in range(1,33,2)]
    phase = [0] + phase + [1]

    red_norm = [ k/255 for k in red ] + [ red[0]/255 ]
    green_norm = [ k/255 for k in green ] + [ green[0]/255 ]
    blue_norm = [ k/255 for k in blue ] + [ blue[0]/255 ]

    redtuple=[]
    greentuple=[]
    bluetuple=[]
    for j in range(18):
        redtuple.append((phase[j],red_norm[j],red_norm[j+1]))
        greentuple.append((phase[j],green_norm[j],green_norm[j+1]))
        bluetuple.append((phase[j],blue_norm[j],blue_norm[j+1]))

    redtuple=tuple(redtuple)
    greentuple=tuple(greentuple)
    bluetuple=tuple(bluetuple)

    cdict = { 'red': redtuple, 'green': greentuple, 'blue': bluetuple }

    return cdict

#%%
def make_geotiff(data, length, width, latn, lonw, dlat, dlon, outfile):
    ## Grid registration to pixel registration by shifing half pixel
    latn_p = latn - dlat/2
    lonw_p = lonw - dlon/2

    driver = gdal.GetDriverByName('GTiff')
    outRaster = driver.Create(outfile, width, length, 1, gdal.GDT_Float32)
    outRaster.SetGeoTransform((lonw_p, dlon, 0, latn_p, 0, dlat))
    outband = outRaster.GetRasterBand(1)
    outband.WriteArray(data)
    outband.SetNoDataValue(0)
    outRaster.SetMetadataItem('AREA_OR_POINT', 'Point')
    outRasterSRS = osr.SpatialReference()
    outRasterSRS.ImportFromEPSG(4326)
    outRaster.SetProjection(outRasterSRS.ExportToWkt())
    outband.FlushCache()

    return


#%% Main
def main(argv=None):
    
    #%% Check argv
    if argv == None:
        argv = sys.argv
        
    start = time.time()
#    print("\n{} {}".format(os.path.basename(argv[0]), ' '.join(argv[1:])), flush=True)


    #%% Set default
    infile = []
    ztdmfile = []
    ztdsfile = []
    LOSufile = []
    dempar = []
    endian = 'little'
    resampleAlg = 'cubicspline'# None # 'cubic' 
    radar_frequency = 5405000000 #Hz, Sentinel-1

    #%% Read options
    try:
        try:
            opts, args = getopt.getopt(argv[1:], "hu:d:m:s:l:p:", ["help", "bigendian"])
        except getopt.error as msg:
            raise Usage(msg)
        for o, a in opts:
            if o == '-h' or o == '--help':
                print(__doc__)
                return 0
            elif o == '-d':
                infile = a
                intype = 'diff'
            elif o == '-u':
                infile = a
                intype = 'unw'
            elif o == '-m':
                ztdmfile = a
                ztdmpar = ztdmfile+'.rsc'
            elif o == '-s':
                ztdsfile = a
            elif o == '-l':
                LOSufile = a
            elif o == '-p':
                dempar = a
            elif o == '--bigendian':
                endian = 'big'

        if not infile:
            raise Usage('No input file given, use either -u or -d!')
        elif not os.path.exists(infile):
            raise Usage('No {} exists!'.format(infile))
        if not ztdmfile:
            raise Usage('No ztd file for master given, -m is not optional!')
        elif not os.path.exists(ztdmfile):
            raise Usage('No {} exists!'.format(ztdmfile))
        elif not os.path.exists(ztdmpar):
            raise Usage('No {} exists!'.format(ztdmpar))
        if not ztdsfile:
            raise Usage('No ztd file for slave given, -s is not optional!')
        elif not os.path.exists(ztdsfile):
            raise Usage('No {} exists!'.format(ztdsfile))
        if os.path.getsize(ztdmfile)!=os.path.getsize(ztdsfile):
            raise Usage('File size of two ztd files are not same!')
        if not LOSufile:
            raise Usage('No incidence angle given, -l is not optional!')
        if dempar and not os.path.exists(dempar):
            raise Usage('No {} exists!'.format(dempar))

    except Usage as err:
        print("\nERROR:", file=sys.stderr, end='')
        print("  "+str(err.msg), file=sys.stderr)
        print("\nFor help, use -h or --help.\n", file=sys.stderr)
        return 2
    
    
    #%% Read data information
    ### Get general info
    speed_of_light = 299792458 #m/s
    wavelength = speed_of_light/radar_frequency #meter
    m2r_coef = 4*np.pi/wavelength

    ### Read infile. Grid registration
    if os.path.splitext(infile)[1] == '.tif': ## GeoTIFF
        informat = 'GeoTIFF'
        geotiff = gdal.Open(infile)
        phase = geotiff.ReadAsArray()
        width_geo = geotiff.RasterXSize
        length_geo = geotiff.RasterYSize
        lonw_p, dlon_geo, _, latn_p, _, dlat_geo = geotiff.GetGeoTransform()
        ## lat lon are in pixel registration. dlat is negative
        latn_geo = latn_p + dlat_geo/2
        lonw_geo = lonw_p + dlon_geo/2
        lats_geo = latn_geo+dlat_geo*(length_geo-1)
        lone_geo = lonw_geo+dlon_geo*(width_geo-1)
    else: ## float32 or complex64
        if not dempar :
            print('\nERROR: Input file seems not to be GeoTIFF, -p option is required!', file=sys.stderr)
            return 1
        
        width_geo = int(get_param_par(dempar, 'width'))
        length_geo = int(get_param_par(dempar, 'nlines'))
        dlat_geo = float(get_param_par(dempar, 'post_lat')) #minus
        dlon_geo = float(get_param_par(dempar, 'post_lon'))
        latn_geo = float(get_param_par(dempar, 'corner_lat'))
        lonw_geo = float(get_param_par(dempar, 'corner_lon'))
        lats_geo = latn_geo+dlat_geo*(length_geo-1)
        lone_geo = lonw_geo+dlon_geo*(width_geo-1)
        
        if os.path.getsize(infile) == width_geo*length_geo*4:
            informat = 'float32'
            phase = read_img(infile, length_geo, width_geo, dtype=np.float32, endian=endian)
        elif os.path.getsize(infile) == width_geo*length_geo*8:
            informat = 'complex64'
            data = read_img(infile, length_geo, width_geo, dtype=np.complex64, endian=endian)
            amp = np.abs(data)
            phase = np.angle(data)
        else:
            print('\nERROR: {} seems to have wrong dimension!'.format(infile), file=sys.stderr)
            return 1
            
    print('\nFile format of {}: {}'.format(infile, informat))
    print('File type   of {}: {}'.format(infile, intype))

    phase[phase==0] = np.nan ## 0 means no data

    ### Read ztd. Grid registration
    width_ztd = int(get_param_par(ztdmpar, 'WIDTH'))
    length_ztd = int(get_param_par(ztdmpar, 'FILE_LENGTH'))
    dlat_ztd = float(get_param_par(ztdmpar, 'Y_STEP')) #minus
    dlon_ztd = float(get_param_par(ztdmpar, 'X_STEP'))
    latn_ztd = float(get_param_par(ztdmpar, 'Y_FIRST'))
    lonw_ztd = float(get_param_par(ztdmpar, 'X_FIRST'))


    ### Read LOSu and calc incidence angle
    if os.path.exists(LOSufile):
        if os.path.splitext(LOSufile)[1] == '.tif': ## GeoTIFF
            print('\nGeoTIFF file of vertical component of LOS vector: {}'.format(LOSufile))
            LOSu = gdal.Open(LOSufile).ReadAsArray()
            if LOSu.shape != (length_geo, width_geo):
                print('\nERROR: File size of {} is not same as {}!'.format(LOSufile, infile), file=sys.stderr)
                return 1
        else: ## float32
            print('\nfloat32 file of vertical component of LOS vector: {}'.format(LOSufile))
            LOSu = read_img(LOSufile, length_geo, width_geo, dtype=np.float32, endian=endian)

        LOSu[LOSu==0] = np.nan

    else:
        try:
            print('\nConstant incidence angle (degree): {}'.format(LOSufile))
            LOSu = np.cos(float(LOSufile)*np.pi/180)
        except:
            print('\nERROR: {} must be a GeoTIFF file or scholar!'.format(LOSufile), file=sys.stderr)
            return 1


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
    hdrmfile = ztdmfile+'.hdr'
    hdrsfile = ztdsfile+'.hdr'
    with open(hdrmfile, "w") as f:
        f.write("\n".join(strings))
    shutil.copyfile(hdrmfile, hdrsfile)


    #%% Process ztd files 
    bilmfile = ztdmfile+'.bil'
    if os.path.exists(bilmfile): os.remove(bilmfile)
    os.symlink(ztdmfile, bilmfile)

    bilsfile = ztdsfile+'.bil'
    if os.path.exists(bilsfile): os.remove(bilsfile)
    os.symlink(ztdsfile, bilsfile)

    ### Cut and resapmle ztd to geo
    print('\nCut and resapmle 2 ztd...', flush=True)
    ztdm_geo = gdal.Warp("", bilmfile, format='MEM', outputBounds=(lonw_geo, lats_geo, lone_geo, latn_geo), width=width_geo, height=length_geo, resampleAlg=resampleAlg, srcNodata=0).ReadAsArray()
    ztds_geo = gdal.Warp("", bilsfile, format='MEM', outputBounds=(lonw_geo, lats_geo, lone_geo, latn_geo), width=width_geo, height=length_geo, resampleAlg=resampleAlg, srcNodata=0).ReadAsArray()

    ### differencem, Meter to rad, ztd->sltd
    print('\nCompute dsltd from 2 ztd...', flush=True)
    dsltd = (ztds_geo-ztdm_geo)*m2r_coef/LOSu ## LOSu=cos(inc)
    dsltd[dsltd==0] = np.nan

    ### Correct phase
    print('\nCorrect phase...', flush=True)
    phase_cor = phase-dsltd
    
    
    if intype == 'diff': ## Wrap
        phase_cor = np.angle(np.exp(1j*phase_cor))
        dsltd = np.angle(np.exp(1j*dsltd))
    else: ## Calc std if unw
        std_phase = np.nanstd(phase[~np.isnan(phase_cor)])
        std_phasecor = np.nanstd(phase_cor)
        rate = (std_phase-std_phasecor)/std_phase*100
        print('\nPhase STD Before Correction: {:4.1f} rad'.format(std_phase))
        print('Phase STD After  Correction: {:4.1f} rad'.format(std_phasecor))
        print('Reduction Rate             : {:5.1f} %'.format(rate))

    os.remove(hdrmfile)
    os.remove(hdrsfile)
    os.remove(bilmfile)
    os.remove(bilsfile)

    
    #%% Output dsltd and phase corrected
    print('\nOutput files...', flush=True)
    if informat == 'GeoTIFF':
        dsltdfile = infile.replace('.tif', '.dsltd.tif')
        outfile = infile.replace('.tif', '.gacos.tif')
        pngfile = infile.replace('.tif', '.gacos.png')

        make_geotiff(dsltd, length_geo, width_geo, latn_geo, lonw_geo, dlat_geo, dlon_geo, dsltdfile)
        make_geotiff(phase_cor, length_geo, width_geo, latn_geo, lonw_geo, dlat_geo, dlon_geo, outfile)
        
    else: ## float32 or complex64
        dsltdfile = infile+'.dsltd'
        outfile = infile+'.gacos'
        pngfile = infile+'.gacos.png'

        write_img(dsltd, dsltdfile, endian)
        if informat == 'complex64':
            data = amp*np.exp(1j*phase_cor)
            data[np.isnan(phase_cor)] = 0
            write_img(data, outfile, endian)
        else: ## float32
            write_img(phase_cor, outfile, endian)
        
    
    #%% Output png for comparison
    print('\nOutput png image...', flush=True)
    if intype == 'diff':
        data3 = [phase, phase_cor, dsltd]
        title3 = ['Phase', 'Corrected phase', 'Correction']
    else: ## unw
        cycle = 3  # 2pi*3/cycle for png
        data3 = [np.angle(np.exp(1j*(data/cycle))*cycle) for data in [phase, phase_cor, dsltd]]
        title3 = ['Phase (STD: {:.1f} rad)'.format(std_phase), 'Corrected phase (STD: {:.1f} rad)'.format(std_phasecor), 'Correction ({:.1f}% reduced)'.format(rate)]
    
    make_3im_png(data3, pngfile, 'insar', title3, vmin=-np.pi, vmax=np.pi, cbar=False)
        
    
    #%% Finish
    elapsed_time = time.time()-start
    hour = int(elapsed_time/3600)
    minite = int(np.mod((elapsed_time/60),60))
    sec = int(np.mod(elapsed_time,60))
    print("\nElapsed time: {0:02}h {1:02}m {2:02}s".format(hour,minite,sec))

    print('\n{} Successfully finished!!\n'.format(os.path.basename(argv[0])))
    print('Output files:')
    print('  {} (corrected data)'.format(outfile))
    print('  {} (comparison figure)'.format(pngfile))
    print('  {} (correction data)\n'.format(dsltdfile))


#%% main
if __name__ == "__main__":
    sys.exit(main())
