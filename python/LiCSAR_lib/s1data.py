#!/usr/bin/env python

import LiCSquery as lq
import os
import datetime as dt
from sentinelsat.sentinel import SentinelAPI
import re
import arch2DB
#here is the nostdout function:
from LiCSAR_misc import *

def get_images_for_frame(frameName, startdate = dt.datetime.strptime('20141001','%Y%m%d').date(), enddate = dt.date.today(), sensType = 'IW'):
    #startdate and enddate should be of type datetime.date (but datetime may also work?)
    if str(type(startdate)).split("'")[1] == 'str':
        startdate = dt.datetime.strptime(startdate,'%Y%m%d').date()
    if str(type(enddate)).split("'")[1] == 'str':
        enddate = dt.datetime.strptime(enddate,'%Y%m%d').date()
    scihub_user, scihub_pass, scihub_url = get_scihub_creds()
    #sensType = 'IW'
    nmax = 100
    footprint = lq.get_wkt_boundaries(frameName)
    ascdesc = frameName[3]
    if ascdesc == 'D': ascdesc = 'DESCENDING'
    if ascdesc == 'A': ascdesc = 'ASCENDING'
    track = int(frameName[0:3])
    #startdate = dt.datetime.strptime(startdate,'%Y%m%d').date()
    scihub = SentinelAPI(scihub_user, scihub_pass, scihub_url)
    result = scihub.query(footprint, date = (startdate.strftime('%Y%m%d'), enddate.strftime('%Y%m%d')), \
                     platformname = 'Sentinel-1', producttype = 'SLC', \
                     relativeorbitnumber = str(track), sensoroperationalmode = sensType, orbitdirection = ascdesc)
    #scihub.to_geodataframe(result)
    if result:
        images = scihub.to_dataframe(result)['title'].values.tolist()
    else:
        images = None
    return images

def get_scihub_creds():
    scihub_url = "https://scihub.copernicus.eu/dhus"
    credentials = os.path.join(os.environ["LiCSAR_configpath"],'scihub_credentials_wget')
    f = open(credentials,'r')
    scihub_user = re.sub('\W+','',f.readline().split('=')[1])
    scihub_pass = re.sub('\W+','',f.readline().split('=')[1])
    f.close()
    return scihub_user, scihub_pass, scihub_url

def get_neodc_path_images(images):
    #should work for both one image or image_list
    imlist = images
    neodc_paths = []
    if str(type(images)).split("'")[1] == 'str':
        imlist = []
        imlist.append(images)
    for image in imlist:
        AorB = image[2].lower()
        year = image[17:21]
        mon = image[21:23]
        day = image[23:25]
        if int(year)>=2019 and int(mon)>=6 and int(day)>=25:
            vers='3'
        else:
            vers='2'
        neodcpath = os.path.join('/neodc/sentinel1'+AorB,'data/IW/L1_SLC/IPF_v'+vers,year,mon,day,image+'.zip')
        neodc_paths.append(neodcpath)
    return neodc_paths

def import_to_licsinfo(images, meta = True):
    #output is list of files to be downloaded
    todown=[]
    neodc_paths = get_neodc_path_images(images)
    for imagepath in neodc_paths:
        if os.path.exists(imagepath):
            arch2DB.main('-f'+imagepath)
        else:
            metaonly = imagepath.replace('.zip','.metadata_only.zip')
            if os.path.exists(metaonly) and meta==True:
                arch2DB.main('-f'+metaonly)
            else:
                print('Image '+os.path.basename(imagepath)+' is not in neodc. Including to the list for download')
                todown.append(os.path.basename(imagepath))
    return todown

def check_and_import_to_licsinfo(frameName, startdate = dt.datetime.strptime('20141001','%Y%m%d').date(), enddate = dt.date.today(), meta = True):
    print('Checking S1 images available on /neodc and importing them to licsinfo database')
    images = get_images_for_frame(frameName, startdate, enddate)
    if images:
        print('There are {0} acquired images within the given period'.format(str(len(images))))
    else:
        print('No images found in scihub for the given dates between {0} and {1}. Aborting'.format(str(startdate),str(enddate)))
        return None
    print('Checking and importing to LiCSInfo database')
    with nostdout():
        todown = import_to_licsinfo(images, meta)
    print('There are '+str(len(todown))+' images acquired that do not exist on /neodc for the given period')
    return todown
