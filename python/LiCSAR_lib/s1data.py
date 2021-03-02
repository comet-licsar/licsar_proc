#!/usr/bin/env python

import LiCSquery as lq
import os
import datetime as dt
from sentinelsat.sentinel import SentinelAPI
import re
#import arch2DB
#here is the nostdout function:
from LiCSAR_misc import *

def download(uuid, slcdir):
    scihub_user, scihub_pass, scihub_url = get_scihub_creds()
    scihub = SentinelAPI(scihub_user, scihub_pass, scihub_url)
    rc = scihub.download(uuid, slcdir)
    return rc

def get_images_for_frame(frameName, startdate = dt.datetime.strptime('20141001','%Y%m%d').date(),
             enddate = dt.date.today(), sensType = 'IW', outAspd = False):
    #startdate and enddate should be of type datetime.date (but datetime may also work)
    # problem is that only full days are selected, no search by time
    if str(type(startdate)).split("'")[1] == 'str':
        startdate = dt.datetime.strptime(startdate,'%Y%m%d').date()
    if str(type(enddate)).split("'")[1] == 'str':
        enddate = dt.datetime.strptime(enddate,'%Y%m%d').date()
    #one extra date needed for scihub:
    enddate = enddate + dt.timedelta(days=1)
    #sensType = 'IW'
    nmax = 100
    footprint = lq.get_wkt_boundaries(frameName)
    ascdesc = frameName[3]
    if ascdesc == 'D': ascdesc = 'DESCENDING'
    if ascdesc == 'A': ascdesc = 'ASCENDING'
    track = int(frameName[0:3])
    #startdate = dt.datetime.strptime(startdate,'%Y%m%d').date()
    result = None
    try:
        scihub_user, scihub_pass, scihub_url = get_scihub_creds()
        scihub = SentinelAPI(scihub_user, scihub_pass, scihub_url)
        result = scihub.query(footprint, date = (startdate.strftime('%Y%m%d'), enddate.strftime('%Y%m%d')), \
                     platformname = 'Sentinel-1', producttype = 'SLC', \
                     relativeorbitnumber = str(track), sensoroperationalmode = sensType, orbitdirection = ascdesc)
    except:
        print('error in scihub search, should try CEDA elastic search (no option for detailed search, not using it now..)')
        #TODO!!
        #from elasticsearch import Elasticsearch
        #query = {
        #    "query": {"match_all": {}}
        #    }
        #es = Elasticsearch(["https://elasticsearch.ceda.ac.uk"])
        #es.search(index="ceda-eo", body=query)
        #result = es.indices.get_mapping(index="ceda-eo")
    #scihub.to_geodataframe(result)
    if result:
        if outAspd:
            images = scihub.to_dataframe(result)
        else:
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
        if image.split('_')[1][0] == 'S':
            sensType='SM'
        else:
            sensType='IW'
        if dt.datetime.strptime(year+mon+day,'%Y%m%d') >= dt.datetime.strptime('20190625','%Y%m%d'):
            vers='3'
        else:
            vers='2'
        neodcpath = os.path.join('/neodc/sentinel1'+AorB,'data',sensType,'L1_SLC/IPF_v'+vers,year,mon,day,image+'.zip')
        neodc_paths.append(neodcpath)
    return neodc_paths

def import_to_licsinfo(images, meta = True):
    #output is list of files to be downloaded
    todown=[]
    print('printing from s1data.py for debug')
    neodc_paths = get_neodc_path_images(images)
    i=0
    leng=len(neodc_paths)
    print('there are {} paths to check'.format(leng))
    for imagepath in neodc_paths:
        i=i+1
        print('[{0}/{1}] checking {2}'.format(str(i),str(leng),imagepath))
        if os.path.exists(imagepath):
            # this was not working, donno why
            #arch2DB.main('-f'+imagepath)
            os.system('arch2DB.py -f {}'.format(imagepath))
        else:
            metaonly = imagepath.replace('.zip','.metadata_only.zip')
            print('SLC not existing in neodc. checking '+metaonly)
            if os.path.exists(metaonly) and meta==True:
                #perhaps some long file name problem, donno..
                #arch2DB.main('-f'+metaonly)
                os.system('arch2DB.py -f {}'.format(metaonly))
            else:
                print('Image '+os.path.basename(imagepath)+' is not in neodc. Including to the list for download')
                todown.append(os.path.basename(imagepath))
    return todown

def check_and_import_to_licsinfo(frameName, startdate = dt.datetime.strptime('20141001','%Y%m%d').date(), enddate = dt.date.today(), meta = True):
    print('Checking S1 images available on /neodc and importing them to licsinfo database')
    if frameName.split('_')[1] == 'SM':
        sensType = 'SM'
    else:
        sensType = 'IW'
    images = get_images_for_frame(frameName, startdate, enddate, sensType)
    if images:
        print('There are {0} acquired images within the given period'.format(str(len(images))))
    else:
        print('No images found in scihub for the given dates between {0} and {1}. Aborting'.format(str(startdate),str(enddate)))
        return None
    print('Checking and importing to LiCSInfo database')
    #with nostdout():
    todown = import_to_licsinfo(images, meta)
    if len(todown)>0:
        print('Summary: There are '+str(len(todown))+' images that were physically acquired but do not exist in CEDA Archive')
    else:
        print('All needed files should exist in CEDA Archive')
    return todown
