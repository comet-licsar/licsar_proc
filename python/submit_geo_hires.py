#!/usr/bin/env python
"""
This will prepare higher res geo files needed to generate 50x50 m geocoded ifgs

Usage:
    just run it from the processing folder (e.g. $BATCH_CACHE_DIR/$frame)

"""

from LiCSAR_lib.coreg_lib import geocode_dem
import datetime as dt
import os, glob

outres = 0.0005

framename = os.getcwd().split('/')[-1]
if not (framename.count('_') == 2 and len(framename) == 17):
    print('please make sure you start this script from FRAME folder!')
    print('best if this is your BATCH_CACHE_DIR')
    exit()

procdir = os.getcwd()
geodir = os.path.join(procdir,'geo_50m')
demdir = os.path.join(procdir,'DEM')
masterdate = glob.glob('geo/????????.hgt')[0].split('/')[1].split('.')[0]
masterslcdir = os.path.join(procdir,'SLC',masterdate)
masterdate = dt.datetime.strptime(masterdate,"%Y%m%d")

#os.rename(geodir, 'geo_orig')
os.mkdir(geodir)
os.system('touch {}/locked'.format(geodir))
geocode_dem(masterslcdir,geodir,demdir,procdir,masterdate,outres)
os.system('rm {}/locked'.format(geodir))
