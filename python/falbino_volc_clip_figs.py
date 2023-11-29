#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Proces LiCSAR data for volcanoes

Create required tif and image files for each volcano and frame
"""

# -- imports

# - std lib imports:
from __future__ import division
from ctypes import c_int
import datetime
from glob import glob
import heapq
from itertools import chain
from multiprocessing import Manager, Pool
import os
from subprocess import call
import sys
import time
import warnings

# - third party imports:
import elevation
import great_circle_calculator.great_circle_calculator as gcc
import matplotlib

matplotlib.use('Agg')
import matplotlib.cbook
import matplotlib.pyplot as plt
from matplotlib_scalebar.scalebar import ScaleBar
from mpl_toolkits.basemap import Basemap
import numpy as np
from numpy import arctan
from numpy import arctan2
from numpy import cos
from numpy import gradient
from numpy import pi
from numpy import sin
from numpy import sqrt
from osgeo import osr, gdal

# -- global variables

# region to process from command argument:
COUNTRY = sys.argv[1]
# check region is valid or give up:
if COUNTRY not in [
    'africa', 'atlantic_island', 'central_america', 'eastern_asia', 'europe',
    'iceland', 'indian_island', 'middle_east', 'north_america',
    'northern_asia', 'pacific_island', 'south_america', 'southeast_asia'
]:
    err_msg = 'Unknown region: {0}\n'.format(COUNTRY)
    sys.stderr.write(err_msg)
    sys.exit()

# main volc-proc directory:
HOMEPATH = '/gws/nopw/j04/nceo_geohazards_vol1/projects/LiCS/volc-proc/'
# list database directory:
LISTPATH = HOMEPATH + 'list_database/'
# database directory:
DATABASEPATH = HOMEPATH + 'database_web/'
# path to output directory for this region:
WORKPATH = HOMEPATH + 'current/' + COUNTRY

# path to frame list file:
LIST_FRAME = LISTPATH + COUNTRY + '_frame.txt'
# path to volcano list file:
LIST_VOLCANO = LISTPATH + COUNTRY + '_volcano.txt'
# path to report / log file:
REPORT_FILE = LISTPATH + COUNTRY + '_report.txt'

# distance around volcano (degrees) to crop data:
DX = 0.5
DY = 0.5
# wavelength constants:
WAVELENGTH = 0.055465
RAD2CM = (WAVELENGTH / (4 * pi)) * 100

# multiprocessing pool size:
POOL_SIZE = 4


# --- classes

class SuppressStream(object):
    """
    Class to suppress an output strem (e.g. stderr).
    Method copied from:
      https://stackoverflow.com/a/57677370
    """

    def __init__(self, stream=sys.stderr):
        self.orig_stream_fileno = stream.fileno()
        self.orig_stream_dup = os.dup(self.orig_stream_fileno)
        self.devnull = open(os.devnull, 'w')

    def __enter__(self):
        os.dup2(self.devnull.fileno(), self.orig_stream_fileno)

    def __exit__(self, my_type, my_value, my_traceback):
        os.close(self.orig_stream_fileno)
        os.dup2(self.orig_stream_dup, self.orig_stream_fileno)
        os.close(self.orig_stream_dup)
        self.devnull.close()


# -- functions

def display_progress(progress_count, progress_total):
    """
    Display percentage completion information
    """
    # current completion percentage:
    percent_complete = (progress_count / progress_total) * 100
    # display progress:
    file_msg = '\rProgress : {0:.02f} %'
    sys.stdout.flush()
    sys.stdout.write(file_msg.format(percent_complete))
    sys.stdout.flush()


def mp_wrapper(mp_options):
    """
    Wrapper for multiprocessing function
    """
    # try to catch KeyboardInterrupt:
    try:
        # get options:
        mp_function = mp_options['mp_function']
        mp_count = mp_options['mp_count']
        mp_progress = mp_options['mp_progress']
        mp_lock = mp_options['mp_lock']
        mp_display_progress = mp_options['mp_display_progress']
        # acquire lock and display progress message:
        if mp_display_progress:
            mp_lock.acquire()
            display_progress(mp_progress.value, mp_count)
            mp_lock.release()
        # run the worker function:
        mp_out = mp_function(mp_options)
        # increment progress count:
        mp_progress.value += 1
        # acquire lock and display progress message:
        if mp_display_progress:
            mp_lock.acquire()
            display_progress(mp_progress.value, mp_count)
            mp_lock.release()
        # return the result:
        return mp_out
    except KeyboardInterrupt:
        pass


def mp_run(mp_function, mp_options, pool_size, mp_display_progress=False):
    """
    Run function via multiprocessing pool with specified options and pool
    size
    """
    # create a multiprocessing manager:
    mp_manager = Manager()
    # add progress counter value to the multiprocessing manager:
    mp_progress = mp_manager.Value(c_int, 0)
    # multiprocessing lock for message display:
    mp_lock = mp_manager.Lock()
    # multiprocessing lock for logging:
    mp_log_lock = mp_manager.Lock()
    # update mp_options to add function and multiprocessing manager
    # information:
    for mp_option in mp_options:
        mp_option['mp_function'] = mp_function
        mp_option['mp_count'] = len(mp_options)
        mp_option['mp_progress'] = mp_progress
        mp_option['mp_lock'] = mp_lock
        mp_option['mp_log_lock'] = mp_log_lock
        mp_option['mp_display_progress'] = mp_display_progress
    # create a multiprocessing pool:
    mp_pool = Pool(pool_size)
    # use map_async to run the function:
    mp_result = mp_pool.map_async(mp_wrapper, mp_options)
    # try to exit cleanly on KeyboardInterrupt:
    try:
        mp_out = mp_result.get()
        mp_pool.close()
        mp_pool.join()
    except KeyboardInterrupt:
        # try to kill mp pool nicely:
        try:
            mp_pool.terminate()
            mp_pool.join()
        except:
            pass
        # exit:
        sys.stdout.write('\n')
        sys.exit()
    # update progress and add a line break after execution:
    if mp_display_progress:
        display_progress(1, 1)
        sys.stdout.write('\n')
        sys.stdout.flush()
    # close the multiprocessing manager:
    mp_manager.shutdown()
    # return any result:
    return mp_out


def log_message(message, log_file=None, init_log=False, mp_lock=None):
    """
    Log message to file and also print
    """
    # if no log file defined use global REPORT_FILE ... :
    if not log_file:
        log_file = REPORT_FILE
    # get current date and time for log:
    log_date = datetime.datetime.now()
    log_date_str = log_date.strftime('%Y-%m-%d %H:%M:%S')
    # format the log message:
    log_msg = '{0}: {1}\n'.format(log_date_str, message)
    # if we have a multiprocessing lock, acquire it:
    if mp_lock:
        mp_lock.acquire()
    # if initialising log file:
    if init_log:
        # open log file in write mode:
        log_handle = open(log_file, 'w')
    else:
        # open log file in append mode:
        log_handle = open(log_file, 'a')
    # write the message to file and std out:
    log_handle.write(log_msg)
    sys.stdout.write(log_msg)
    # close the log file:
    log_handle.close()
    # if we have a multiprocessing lock, release it:
    if mp_lock:
        mp_lock.release()


def hillshade(arr, azimuth, angle_altitude):
    """
    Return hillshaded data
    """
    x, y = gradient(arr)
    slope = pi / 2. - arctan(sqrt(x * x + y * y))
    aspect = arctan2(-x, y)
    azimuthrad = azimuth * pi / 180.
    altituderad = angle_altitude * pi / 180.
    shaded = (
            sin(altituderad) * sin(slope) + cos(altituderad) * cos(slope) *
            cos(azimuthrad - aspect)
    )
    return 255 * (shaded + 1) / 2


def get_frame_volcanoes(options):
    """
    Get volcanoes for this frame
    """
    # get options for this process:
    frame = options['frame']
    number = options['number']
    licspath = options['licspath']
    framepath = options['framepath']
    epochpath = options['epochpath']
    paths = options['paths']
    log_lock = options['mp_log_lock']

    # define a log message wrapper which includes lock file:
    def _log_message(message):
        log_message(message, mp_lock=log_lock)

    # init list for storing volcanoes to process:
    frame_out = []

    # count number of diff_pha tif files:
    geotif_files = [
        y for x in os.walk(framepath)
        for y in glob(os.path.join(x[0], '*diff_pha.tif'))
    ]
    nb_lics = (len(geotif_files))
    _log_message('{0} :: Found {1} ifgs'.format(frame, nb_lics))

    # if no tif files, exit:
    if not geotif_files:
        _log_message('{0} :: Frame exists but empty'.format(frame))
        return frame_out

    # get coordinate and projection information for this frame:
    filename = geotif_files[0]
    src = gdal.Open(filename)
    Xmin, xres, xskew, Ymax, yskew, yres = src.GetGeoTransform()
    Xmax = Xmin + (src.RasterXSize * xres)
    Ymin = Ymax + (src.RasterYSize * yres)
    projection = src.GetProjection()
    geotransform = src.GetGeoTransform()
    src = None
    del src

    # init per volcano processing options:
    volc_options = options.copy()
    volc_options['projection'] = projection
    volc_options['geotransform'] = geotransform

    # open volcanoes list for reading:
    volcanoes = open(LIST_VOLCANO, "r")
    # loop through list of volcanoes:
    for columns in (raw.strip().split() for raw in volcanoes):
        volcano = (columns[0].lower())
        xv = float(columns[2])
        yv = float(columns[1])
        # 1- SELECT THE VOLCANO BELONGING TO THE FRAME
        if (xv > (Xmin + DX / 2) and xv < (Xmax - DX / 2) and
                yv > (Ymin + DX / 2) and yv < (Ymax - DX / 2)):
            # log a message about this volcano:
            _log_message('{0} :: Contains {1}'.format(frame, volcano))
            # store options for this volcano:
            my_options = volc_options.copy()
            my_options['volcano'] = volcano
            my_options['xv'] = xv
            my_options['yv'] = yv
            frame_out.append(my_options)
    # close volcanoes file:
    volcanoes.close()
    # log volcano count for frame:
    _log_message('{0} :: Contains {1} volcanoes'.format(frame, len(frame_out)))
    # return the options for all volcanoes in this frame:
    return frame_out


def process_frame_volcano(options):
    """
    Process data for volcano in frame
    """
    # ignore these warnings:
    warnings.filterwarnings('ignore', category=matplotlib.cbook.mplDeprecation)
    # register gdal drivers:
    gdal.AllRegister()

    # get options for this process:
    frame = options['frame']
    number = options['number']
    licspath = options['licspath']
    framepath = options['framepath']
    epochpath = options['epochpath']
    paths = options['paths']
    log_lock = options['mp_log_lock']
    projection = options['projection']
    geotransform = options['geotransform']
    volcano = options['volcano']
    xv = options['xv']
    yv = options['yv']

    # label for this volcano / frame:
    volc_label = '{0}_{1}'.format(volcano, frame)

    # define a log message wrapper which includes lock file:
    def _log_message(message):
        log_message(message, mp_lock=log_lock)

    # init results dict:
    volc_out = {
        'tif': 0,
        'png': 0,
        'jpg': 0
    }

    # set variables for this volcano / frame:
    dirname = volc_label
    volc_dir = os.path.join(WORKPATH, dirname)
    dirtif = os.path.join(WORKPATH, dirname, 'tif')
    dirtif_empty = os.path.join(WORKPATH, dirname, 'empty_files')
    direpoch = os.path.join(WORKPATH, dirname, 'epoch')
    dirpng = os.path.join(WORKPATH, dirname, 'png')
    dirdem = os.path.join(WORKPATH, dirname, 'dem')
    dirmask = os.path.join(WORKPATH, dirname, 'mask')
    # create output directories if required:
    if not os.path.exists(volc_dir):
        os.mkdir(volc_dir)
    if not os.path.exists(dirtif):
        os.mkdir(dirtif)
    if not os.path.exists(dirpng):
        os.mkdir(dirpng)
    if not os.path.exists(dirdem):
        os.mkdir(dirdem)
    if not os.path.exists(dirmask):
        os.mkdir(dirmask)
    if not os.path.exists(direpoch):
        os.mkdir(direpoch)

    # min and max x and y coordinates for this volcano / frame:
    xmin = (xv - DX / 2)
    xmax = (xv + DX / 2)
    ymin = (yv - DY / 2)
    ymax = (yv + DY / 2)

    # 2a- CREATE THE DEM FILE AND HILL SHADE OF THE CROPPING
    demfile1 = dirname + '_SRTM1_dem.tif'
    demfile2 = dirname + '_SRTM1_resize_dem.tif'
    input_dem = os.path.join(dirdem, demfile1)
    output_dem = os.path.join(dirdem, demfile2)
    framepath = licspath + number + '/' + frame + '/metadata/'
    demfile_hgt = framepath + frame + '.geo.hgt.tif'
    if not os.path.isfile(input_dem):
        if ymax > 60 or ymin < -60:
            call([
                'gdalwarp', '-q', '-t_srs', 'EPSG:4326', '-te', str(xmin),
                str(ymin), str(xmax), str(ymax), demfile_hgt, input_dem
            ])
        else:
            # suppress std out when using elevation:
            with SuppressStream(sys.stdout):
                elevation.clip(bounds=(xmin, ymin, xmax, ymax), output=input_dem)
        _log_message('{0} :: DEM ORIGINAL {1} CREATED'.format(volc_label, demfile1))
        volc_out['tif'] += 1
    if not os.path.isfile(output_dem):
        call([
            'gdalwarp', '-q', '-t_srs', 'EPSG:4326', '-te', str(xmin),
            str(ymin), str(xmax), str(ymax), '-ts', str(DX * 1000),
            str(DY * 1000), input_dem, output_dem
        ])
        _log_message('{0} :: DEM RESIZE {1} CREATED'.format(volc_label, demfile2))
        volc_out['tif'] += 1

    # 2b- CREATE MASK BASED ON SENTINEL-2
    dirmask_S2 = '/gws/nopw/j04/nceo_geohazards_vol1/projects/LiCS/proc/current/GIS/'
    maskfile_S2 = 'water_mask.tif'
    input_mask = os.path.join(dirmask_S2, maskfile_S2)
    maskfile2 = dirname + '_mask_water.tif'
    output_mask2 = os.path.join(dirmask, maskfile2)
    if not os.path.isfile(output_mask2):
        call([
            'gdalwarp', '-q', '-t_srs', 'EPSG:4326', '-te', str(xmin),
            str(ymin), str(xmax), str(ymax), '-ts', str(DX * 1000),
            str(DY * 1000), input_mask, output_mask2
        ])
        _log_message('{0} :: MASK {1} CREATED'.format(volc_label, maskfile2))
        volc_out['tif'] += 1

    # 2b- CREATE MASK BASED ON DEM FILE
    maskfile = dirname + '_mask_sea.tif'
    output_mask = os.path.join(dirmask, maskfile)
    if not os.path.isfile(output_mask):
        ds = gdal.Open(output_dem)
        band = ds.GetRasterBand(1)
        rows = ds.RasterYSize
        cols = ds.RasterXSize
        datatype = band.DataType
        arr = band.ReadAsArray()
        driver = ds.GetDriver()
        outds = driver.Create(output_mask, rows, cols, 1, gdal.GDT_Int32)
        outBand = outds.GetRasterBand(1)
        outData = np.ones((rows, cols), np.int16)
        if COUNTRY == 'africa':
            outData[arr < -1000] = 0
        else:
            outData[arr < -0.1] = 0
        outBand.WriteArray(outData, 0, 0)
        outBand.FlushCache()  ##saves to disk!
        outds.SetGeoTransform(ds.GetGeoTransform())
        outds.SetProjection(ds.GetProjection())
        outds = None
        del outData
        ds = None
        del ds
        _log_message('{0} :: MASK {1} CREATED'.format(volc_label, maskfile))
        volc_out['tif'] += 1

    # 3- CROPPING THE TIFF FILE
    for root, dirs, files in chain.from_iterable(
            os.walk(path) for path in paths
    ):
        for file in files:
            if (file.endswith('diff_pha.tif') | file.endswith('cc.tif') |
                    file.endswith('unw.tif') | file.endswith('mli.tif') |
                    file.endswith('sltd.geo.tif') |
                    file.endswith('ztd.geo.tif')):
                filecrop = volcano + '_' + file
                masterdate = filecrop[len(volcano) + 1:len(volcano) + 9]
                slavedate = filecrop[len(volcano) + 10:len(volcano) + 18]
                if (file.endswith('mli.tif') | file.endswith('sltd.geo.tif') |
                        file.endswith('ztd.geo.tif')):
                    output_file = os.path.join(direpoch, filecrop)
                else:
                    output_file = os.path.join(dirtif, filecrop)
                    output_file_empty = os.path.join(dirtif_empty, filecrop)
                filetmp = 'temp_' + volcano + '_' + file
                temp_file = os.path.join(dirtif, filetmp)
                if (not os.path.isfile(output_file) |
                        os.path.isfile(output_file_empty)):
                    input_file = os.path.join(root, file)
                    try:
                        test = gdal.OpenShared(input_file, gdal.GA_ReadOnly)
                    except IOError:
                        _log_message(
                            '{0} :: The File {1} is corrupted'.format(
                                volc_label, file
                            )
                        )
                        continue
                    if test == None:
                        _log_message(
                            '{0} :: The File {1} is corrupted'.format(
                                volc_label, file
                            )
                        )
                        continue
                    prj = test.GetProjection()
                    srs = osr.SpatialReference(wkt=prj)
                    test2 = srs.GetAttrValue('geogcs')
                    if test2 == None:
                        _log_message(
                            '{0} :: The File {1} has an empty LOCAL_CS'.format(
                                volc_label, file
                            )
                        )
                        test.SetGeoTransform(geotransform)
                        test.SetProjection(projection)
                        prj = test.GetProjection()
                        srs = osr.SpatialReference(wkt=prj)
                    data = test.ReadAsArray()
                    data = np.array(data, dtype=np.float64)
                    max_data = np.nanmax(data)
                    min_data = np.nanmin(data)
                    if max_data - min_data == 0:
                        _log_message(
                            '{0} :: The File {1} is empty'.format(
                                volc_label, file
                            )
                        )
                        continue
                    filecrop = volcano + '_' + file
                    masterdate = filecrop[len(volcano) + 1:len(volcano) + 9]
                    slavedate = filecrop[len(volcano) + 10:len(volcano) + 18]
                    if file.endswith('cc.tif'):
                        if max_data > 2:
                            call([
                                'gdalwarp', '-q', '-t_srs', 'EPSG:4326', '-te',
                                str(xmin), str(ymin), str(xmax), str(ymax),
                                '-ts', str(DX * 1000), str(DY * 1000),
                                input_file, temp_file
                            ])
                            call([
                                'gdal_translate', '-q', '-ot', 'float32',
                                '-scale', str(1), str(255), str(0), str(1),
                                temp_file, output_file
                            ])
                            os.remove(temp_file)
                        else:
                            call([
                                'gdalwarp', '-q', '-t_srs', 'EPSG:4326', '-te',
                                str(xmin), str(ymin), str(xmax), str(ymax),
                                '-ts', str(DX * 1000), str(DY * 1000), input_file,
                                output_file
                            ])
                    else:
                        call([
                            'gdalwarp', '-q', '-t_srs', 'EPSG:4326', '-te',
                            str(xmin), str(ymin), str(xmax), str(ymax), '-ts',
                            str(DX * 1000), str(DY * 1000), input_file,
                            output_file
                        ])
                    _log_message(
                        '{0} :: CROPPING {1} DONE'.format(volc_label, filecrop)
                    )
                    volc_out['tif'] += 1
                    # close the input file:
                    test = None
                    del test

    # 3bis- CREATING GACOS TIFF FOR EACH INTERFEROGRAM
    pha_files = [f for f in os.listdir(dirtif) if f.endswith('unw.tif')]
    for file in pha_files:
        GACOS_file = file[:-7] + 'GACOS.tif'
        unwcor_file = file[:-7] + 'unw_GACOS.tif'
        input_unw = os.path.join(dirtif, file)
        output_unw = os.path.join(dirtif, unwcor_file)
        output_GACOS = os.path.join(dirtif, GACOS_file)
        masterdate = file[len(volcano) + 1:len(volcano) + 9]
        slavedate = file[len(volcano) + 10:len(volcano) + 18]
        file1 = volcano + '_' + masterdate + '.sltd.geo.tif'
        file2 = volcano + '_' + slavedate + '.sltd.geo.tif'
        input_file1 = os.path.join(direpoch, file1)
        input_file2 = os.path.join(direpoch, file2)
        if not os.path.isfile(output_GACOS):
            if os.path.isfile(input_file1) & os.path.isfile(input_file2):
                gdal_calc1 = 'gdal_calc.py --quiet -A {0} -B {1} --outfile={2}'
                gdal_calc1 += ' --calc="(B-A)" --NoDataValue=0'
                gdal_calc2 = 'gdal_calc.py --quiet -A {0} -B {1} --outfile={2}'
                gdal_calc2 += ' --calc="(B-A)" --NoDataValue=0'
                gdal_calc1_process = gdal_calc1.format(
                    input_file1, input_file2, output_GACOS
                )
                gdal_calc2_process = gdal_calc2.format(
                    output_GACOS, input_unw, output_unw
                )
                os.system(gdal_calc1_process)
                _log_message(
                    '{0} :: CREATING {1} DONE'.format(volc_label, GACOS_file)
                )
                volc_out['tif'] += 1
                os.system(gdal_calc2_process)
                _log_message(
                    '{0} :: CREATING {1} DONE'.format(volc_label, unwcor_file)
                )
                volc_out['tif'] += 1

    # 4- CREATE THE PNG FILE
    pha_files = [f for f in os.listdir(dirtif) if f.endswith('unw.tif')]
    for file in pha_files:
        pngfile = file[:-7] + 'pha_unw_coh.png'
        output_png = os.path.join(dirpng, pngfile)
        masterdate = pngfile[len(volcano) + 1:len(volcano) + 9]
        slavedate = pngfile[len(volcano) + 10:len(volcano) + 18]
        titlename = (
                masterdate[0:4] + '-' + masterdate[4:6] + '-' + masterdate[6:8] +
                '    ' + slavedate[0:4] + '-' + slavedate[4:6] + '-' + slavedate[6:8]
        )
        pngfile2 = file[:-7] + 'GACOS_correction.png'
        output_png2 = os.path.join(dirpng, pngfile2)

        if not os.path.isfile(output_png):
            ds = gdal.Open(output_dem)
            band = ds.GetRasterBand(1)
            array = band.ReadAsArray()
            dem_array = hillshade(array, 315, 45)
            if COUNTRY == 'africa':
                dem_array[array < -1000] = np.nan
            else:
                dem_array[array < 0] = np.nan
            dataset = gdal.Open(output_file)
            data = dataset.ReadAsArray()
            proj_wkt = dataset.GetProjection()
            proj = osr.SpatialReference()
            proj.ImportFromWkt(proj_wkt)
            extent_lonlat = (
                xmin,
                xmax,
                ymin,
                ymax
            )
            # plot the data
            map_plot = Basemap(projection='cyl', resolution='l',
                               llcrnrlon=extent_lonlat[0], llcrnrlat=extent_lonlat[2],
                               urcrnrlon=extent_lonlat[1], urcrnrlat=extent_lonlat[3])
            # plot dem
            val_dpi = 200
            ptA = (xv, ymin)
            ptB = (xv + 1, ymin)
            xscale = gcc.distance_between_points(ptA, ptB)
            #
            coh_file = file[:-7] + 'cc.tif'
            pha_file = file[:-7] + 'diff_pha.tif'
            unw_file = file[:-7] + 'unw.tif'
            path_file = os.path.join(dirtif, coh_file)
            #
            if not os.path.isfile(path_file):
                continue
            test = gdal.Open(path_file)
            coh_data = test.ReadAsArray()
            coh_data = np.array(coh_data, dtype=np.float64)
            path_file = os.path.join(dirtif, pha_file)
            if not os.path.isfile(path_file):
                continue
            test = gdal.Open(path_file)
            pha_data = test.ReadAsArray()
            pha_data = np.array(pha_data, dtype=np.float64)
            path_file = os.path.join(dirtif, unw_file)
            if not os.path.isfile(path_file):
                continue
            test = gdal.Open(path_file)
            unw_data = test.ReadAsArray()
            unw_data = np.array(unw_data, dtype=np.float64)
            #
            parallels = np.arange(-180, 180, 0.25)
            meridians = np.arange(0, 360, 0.25)
            suppress_ticks = False
            #
            fig = plt.figure()
            plt.rcParams.update({'font.size': 6})
            if COUNTRY == 'africa':
                coh_data[array < -1000] = np.nan
                pha_data[array < -1000] = np.nan
                unw_data[array < -1000] = np.nan
            else:
                coh_data[array < 0] = np.nan
                pha_data[array < 0] = np.nan
                unw_data[array < 0] = np.nan
            #
            pha_data = pha_data * RAD2CM
            pha_data[pha_data == 0.0] = np.nan
            #
            ax1 = fig.add_subplot(131)
            ax1.set_title('Wrapped LOS change [cm]', size=8)
            map_plot.drawparallels(
                parallels, labels=[1, 0, 0, 0], linewidth=0.5, dashes=[6, 900]
            )
            map_plot.drawmeridians(
                meridians, labels=[0, 0, 0, 1], linewidth=0.5, dashes=[6, 900]
            )
            suppress_ticks = False
            im1_dem = plt.imshow(
                dem_array, cmap='Greys', alpha=1, extent=extent_lonlat
            )
            im1 = plt.imshow(
                pha_data, cmap='RdYlBu_r', vmin=-1.4, vmax=1.4, alpha=0.5,
                extent=extent_lonlat
            )
            cbar = map_plot.colorbar(
                im1, ax=ax1, location='bottom', pad=0.2, alpha=1
            )
            cbar.solids.set(alpha=1)
            scalebar = ScaleBar(
                xscale, units="m", location="lower right",
                length_fraction=0.25, sep=1, border_pad=0.2
            )
            plt.gca().add_artist(scalebar)
            #
            unw_data = unw_data * RAD2CM
            unw_data[unw_data == 0.0] = np.nan
            unw_mean = np.nanmean(unw_data)
            unw_data = unw_data - unw_mean
            unw_min = np.nanmin(unw_data)
            unw_max = np.nanmax(unw_data)
            unw_std = np.nanstd(unw_data)
            v_unw = max(abs(unw_min), abs(unw_max))
            ax2 = fig.add_subplot(132)
            ax2.set_title('Unwrapped LOS change [cm]', size=8)
            #
            im2_dem = plt.imshow(
                dem_array, cmap='Greys', alpha=1, extent=extent_lonlat
            )
            im2 = plt.imshow(
                unw_data, cmap='RdYlBu_r', vmin=-2 * unw_std, vmax=2 * unw_std,
                alpha=0.5, extent=extent_lonlat
            )
            #
            map_plot.drawparallels(
                parallels, linewidth=0.5, dashes=[6, 900]
            )
            map_plot.drawmeridians(
                meridians, linewidth=0.5, dashes=[6, 900]
            )
            #
            ax2.tick_params(labelbottom=False)
            cbar = map_plot.colorbar(
                im2, ax=ax2, alpha=1, location="bottom", pad=0.2
            )
            cbar.solids.set(alpha=1)
            #
            coh_data[coh_data == 0.0] = np.nan
            ax3 = fig.add_subplot(133)
            ax3.set_title('Coherence', size=8)
            im3_dem = plt.imshow(
                dem_array, cmap='Greys', alpha=1, extent=extent_lonlat
            )
            im3 = plt.imshow(
                coh_data, cmap='viridis', vmin=0, vmax=1, alpha=0.5,
                extent=extent_lonlat
            )
            #
            map_plot.drawparallels(parallels, linewidth=0.5, dashes=[6, 900])
            map_plot.drawmeridians(meridians, linewidth=0.5, dashes=[6, 900])
            ax3.tick_params(labelbottom=False)
            ax3.tick_params(labelleft=False)
            cbar = map_plot.colorbar(
                im3, ax=ax3, alpha=1, location="bottom", pad=0.2
            )
            cbar.solids.set(alpha=1)
            #
            fig = ax3.get_figure()
            fig.tight_layout()
            fig.subplots_adjust(top=0.95)
            fig.suptitle(titlename, y=0.8, size=10)
            # save file
            plt.savefig(
                output_png, format='png', dpi=val_dpi, bbox_inches='tight',
                pad_inches=0.02
            )
            plt.clf()
            plt.cla()
            plt.close()
            ds = None
            dataset = None
            test = None
            del ds, dataset, test
            _log_message('{0} :: PNG {1} CREATED'.format(volc_label, pngfile))
            volc_out['png'] += 1

        if not os.path.isfile(output_png2):
            ds = gdal.Open(output_dem)
            band = ds.GetRasterBand(1)
            array = band.ReadAsArray()
            dem_array = hillshade(array, 315, 45)
            if COUNTRY == 'africa':
                dem_array[array < -1000] = np.nan
            else:
                dem_array[array < 0] = np.nan
            dataset = gdal.Open(output_file)
            data = dataset.ReadAsArray()
            proj_wkt = dataset.GetProjection()
            proj = osr.SpatialReference()
            proj.ImportFromWkt(proj_wkt)
            extent_lonlat = (
                xmin,
                xmax,
                ymin,
                ymax
            )
            # plot the data
            map_plot = Basemap(projection='cyl', resolution='l',
                               llcrnrlon=extent_lonlat[0], llcrnrlat=extent_lonlat[2],
                               urcrnrlon=extent_lonlat[1], urcrnrlat=extent_lonlat[3])
            val_dpi = 200
            ptA = (xv, ymin)
            ptB = (xv + 1, ymin)
            xscale = gcc.distance_between_points(ptA, ptB)
            #
            unw_file = file[:-7] + 'unw.tif'
            GACOS_file = file[:-7] + 'GACOS.tif'
            unwcor_file = file[:-7] + 'unw_GACOS.tif'
            #
            path_file = os.path.join(dirtif, unwcor_file)
            if not os.path.isfile(path_file):
                continue
            test = gdal.Open(path_file)
            unwcor_data = test.ReadAsArray()
            unwcor_data = np.array(unwcor_data, dtype=np.float64)
            #
            path_file = os.path.join(dirtif, unw_file)
            if not os.path.isfile(path_file):
                continue
            test = gdal.Open(path_file)
            unw_data = test.ReadAsArray()
            unw_data = np.array(unw_data, dtype=np.float64)
            #
            path_file = os.path.join(dirtif, GACOS_file)
            if not os.path.isfile(path_file):
                continue
            test = gdal.Open(path_file)
            GACOS_data = test.ReadAsArray()
            GACOS_data = np.array(GACOS_data, dtype=np.float64)
            parallels = np.arange(-180, 180, 0.25)
            meridians = np.arange(0, 360, 0.25)
            suppress_ticks = False
            #
            fig = plt.figure()
            plt.rcParams.update({'font.size': 6})
            if COUNTRY == 'africa':
                unwcor_data[array < -1000] = np.nan
                unw_data[array < -1000] = np.nan
                GACOS_data[array < -1000] = np.nan
            else:
                unwcor_data[array < 0] = np.nan
                unw_data[array < 0] = np.nan
                GACOS_data[array < 0] = np.nan
            #
            unw_data = unw_data * RAD2CM
            unw_data[unw_data == 0.0] = np.nan
            unw_mean = np.nanmean(unw_data)
            unw_data = unw_data - unw_mean
            unw_min = np.nanmin(unw_data)
            unw_max = np.nanmax(unw_data)
            unw_std = np.nanstd(unw_data)
            v_unw = max(abs(unw_min), abs(unw_max))
            #
            ax1 = fig.add_subplot(131)
            ax1.set_title('Unwrapped LOS change [cm]', size=8)
            map_plot.drawparallels(
                parallels, labels=[1, 0, 0, 0], linewidth=0.5, dashes=[6, 900]
            )
            map_plot.drawmeridians(
                meridians, labels=[0, 0, 0, 1], linewidth=0.5, dashes=[6, 900]
            )
            suppress_ticks = False
            im1_dem = plt.imshow(
                dem_array, cmap='Greys', alpha=1, extent=extent_lonlat
            )
            im1 = plt.imshow(
                unw_data, cmap='RdYlBu_r', vmin=-2 * unw_std, vmax=2 * unw_std,
                alpha=0.5, extent=extent_lonlat
            )
            cbar = map_plot.colorbar(
                im1, ax=ax1, location='bottom', pad=0.2, alpha=1
            )
            cbar.solids.set(alpha=1)
            scalebar = ScaleBar(
                xscale, units="m", location="lower right",
                length_fraction=0.25, sep=1, border_pad=0.2
            )
            plt.gca().add_artist(scalebar)
            #
            GACOS_data = GACOS_data * RAD2CM
            GACOS_data[GACOS_data == 0.0] = np.nan
            GACOS_mean = np.nanmean(GACOS_data)
            GACOS_data = GACOS_data - GACOS_mean
            GACOS_min = np.nanmin(GACOS_data)
            GACOS_max = np.nanmax(GACOS_data)
            v_GACOS = max(abs(GACOS_min), abs(GACOS_max))
            #
            ax2 = fig.add_subplot(132)
            ax2.set_title('GACOS LOS change [cm]', size=8)
            im2_dem = plt.imshow(
                dem_array, cmap='Greys', alpha=1, extent=extent_lonlat
            )
            im2 = plt.imshow(
                GACOS_data, cmap='RdYlBu_r', vmin=-2 * unw_std, vmax=2 * unw_std,
                alpha=0.5, extent=extent_lonlat
            )
            map_plot.drawparallels(parallels, linewidth=0.5, dashes=[6, 900])
            map_plot.drawmeridians(meridians, linewidth=0.5, dashes=[6, 900])
            ax2.tick_params(labelbottom=False)
            cbar = map_plot.colorbar(
                im2, ax=ax2, alpha=1, location="bottom", pad=0.2
            )
            cbar.solids.set(alpha=1)
            #
            unwcor_data = unwcor_data * RAD2CM
            unwcor_data[unwcor_data == 0.0] = np.nan
            unwcor_mean = np.nanmean(unwcor_data)
            unwcor_data = unwcor_data - unwcor_mean
            unwcor_min = np.nanmin(unwcor_data)
            unwcor_max = np.nanmax(unwcor_data)
            v_unwcor = max(abs(unwcor_min), abs(unwcor_max))
            #
            ax3 = fig.add_subplot(133)
            ax3.set_title('Corrected Unwrapped LOS change [cm]', size=8)
            im3_dem = plt.imshow(
                dem_array, cmap='Greys', alpha=1, extent=extent_lonlat
            )
            im3 = plt.imshow(
                unwcor_data, cmap='RdYlBu_r', vmin=-2 * unw_std, vmax=2 * unw_std,
                alpha=0.5, extent=extent_lonlat
            )
            map_plot.drawparallels(parallels, linewidth=0.5, dashes=[6, 900])
            map_plot.drawmeridians(meridians, linewidth=0.5, dashes=[6, 900])
            ax3.tick_params(labelbottom=False)
            ax3.tick_params(labelleft=False)
            cbar = map_plot.colorbar(
                im3, ax=ax3, alpha=1, location="bottom", pad=0.2
            )
            cbar.solids.set(alpha=1)
            #
            fig = ax3.get_figure()
            fig.tight_layout()
            fig.subplots_adjust(top=0.95)
            fig.suptitle(titlename, y=0.8, size=10)
            #
            plt.savefig(
                output_png2, format='png', dpi=val_dpi, bbox_inches='tight',
                pad_inches=0.02
            )
            plt.clf()
            plt.cla()
            plt.close()
            ds = None
            dataset = None
            test = None
            del ds, dataset, test
            _log_message(
                '{0} :: PNG GACOS {1} CREATED'.format(volc_label, pngfile2)
            )
            volc_out['png'] += 1

    # 4- CREATE JPG FILE AND TRANSFERT TO DEFOWEB
    if not volc_out['png'] == 0:
        volcano2 = volcano.replace('_', '-')
        defoweb = DATABASEPATH + 'figures_sentinel/' + volcano + '/'
        if not os.path.exists(defoweb):
            os.makedirs(defoweb)
        today = datetime.datetime.today()
        pha_files = [
            f for f in os.listdir(dirpng) if f.endswith('pha_unw_coh.png')
        ]
        masterdate = [
            x[len(volcano) + 1:(len(volcano) + 9)] for x in pha_files
        ]
        slavedate = [
            y[(len(volcano) + 10):(len(volcano) + 18)] for y in pha_files
        ]
        diff1 = [
            (today - (datetime.datetime.strptime(i, '%Y%m%d'))).days
            for i in masterdate
        ]
        diff2 = [
            (today - (datetime.datetime.strptime(i, '%Y%m%d'))).days
            for i in slavedate
        ]
        diff = [diff1[i] + diff2[i] for i in range(len(diff1))]
        nlesser = heapq.nsmallest(3, diff)
        ind = np.argsort(diff)
        Bt = [diff1[i] - diff2[i] for i in range(len(diff1))]
        if len(nlesser) == 3:
            file1_pha = os.path.join(dirpng, pha_files[ind[0]])
            file2_pha = os.path.join(dirpng, pha_files[ind[1]])
            file3_pha = os.path.join(dirpng, pha_files[ind[2]])
            jpgfile = volcano2 + '_' + frame + '.jpg'
            output_jpg = os.path.join(defoweb, jpgfile)
            montage_cmd = 'montage {0} {1} {2} -tile 1x3 -geometry +0+0'
            montage_cmd += ' -resize %40 {3}'
            montage_cmd = montage_cmd.format(
                file1_pha, file2_pha, file3_pha, output_jpg
            )
            os.system(montage_cmd)
            _log_message(
                '{0} :: JPG {1} CREATED'.format(volc_label, jpgfile)
            )
            volc_out['jpg'] += 1

    # log file creation counts:
    _log_message(
        '{0} :: {1} GEOTIFF, {2} PNG and {3} JPG have been created'.format(
            volc_label, volc_out['tif'], volc_out['png'], volc_out['jpg']
        ))
    # return the results:
    return volc_out


# --

def main():
    """
    Main program function
    """
    # init log and log a start up message:
    log_message('Starting LiCSAR-Volcano on {0}'.format(COUNTRY), init_log=True)
    # get a starting time:
    start_time = time.time()
    # init list for multiprocessing options:
    mp_options = []
    # open frames file for reading:
    frames = open(LIST_FRAME, "r")

    # loop through frames:
    for columns in (raw.strip().split() for raw in frames):
        # frame name:
        frame = columns[0]
        # get frame path number:
        number = int(frame[:3])
        number = str(number)

        # OPTION 1 - download on public
        licspath = '/gws/nopw/j04/nceo_geohazards_vol1/public/LiCSAR_products/'
        framepath = licspath + number + '/' + frame + '/' + 'interferograms/'
        epochpath = licspath + number + '/' + frame + '/' + 'epochs/'
        paths = (framepath, epochpath)

        # check if frame exists:
        if os.path.exists(framepath):
            log_message('{0} :: Frame exists'.format(frame))
            # create options for this frame and store:
            mp_options.append({
                'frame': frame,
                'number': number,
                'licspath': licspath,
                'framepath': framepath,
                'epochpath': epochpath,
                'paths': paths
            })
        # move on if frame directory does not exist:
        else:
            log_message('{0} :: Frame does not exist'.format(frame))

    # close frames file:
    frames.close()

    # how many frames we have found:
    frame_count = len(mp_options)

    # get volcanoes in each frame using multiprocessing:
    log_message('Searching for volcanoes in {0} frames'.format(frame_count))
    frames_volcanoes = mp_run(get_frame_volcanoes, mp_options, POOL_SIZE)
    # join multiprocessing results in to single list:
    frames_volcanoes = sum(frames_volcanoes, [])
    # count of volcanoes to process:
    volcano_count = len(frames_volcanoes)
    log_message('Found {0} volcanoes in {1} frames'.format(
        volcano_count, frame_count
    ))

    # process volcanoes using multiprocessing:
    volcanoes_out = mp_run(process_frame_volcano, frames_volcanoes, POOL_SIZE)
    # get file creation counts:
    tif_count = sum([i['tif'] for i in volcanoes_out])
    png_count = sum([i['png'] for i in volcanoes_out])
    jpg_count = sum([i['jpg'] for i in volcanoes_out])
    log_message('{0} GEOTIFF, {1} PNG and {2} JPG have been created'.format(
        tif_count, png_count, jpg_count
    ))

    # clean up elevation data:
    elevation.clean()

    # log completion message:
    run_time = round(time.time() - start_time)
    log_message(
        'Finishing LiCSAR-Volcano on {0} in {1} seconds'.format(
            COUNTRY, run_time
        )
    )
    # all done.


if __name__ == '__main__':
    # try to catch KeyboardInterrupt:
    try:
        main()
    except KeyboardInterrupt:
        sys.stdout.write('\n')
        sys.exit()
