#!/usr/bin/env python

import os
import datetime as dt
from sentinelsat.sentinel import SentinelAPI
import re
import pandas as pd
import requests
#import arch2DB
#here is the nostdout function:
from LiCSAR_misc import *

def get_info_pd(fileid = 'S1A_IW_SLC__1SDV_20210908T235238_20210908T235305_039597_04AE3C_4CA7', returncol = None):
    filename = fileid+'.SAFE'
    scihub_user, scihub_pass, scihub_url = get_scihub_creds()
    scihub = SentinelAPI(scihub_user, scihub_pass, scihub_url)
    info = scihub.query(filename = filename)
    info = scihub.to_dataframe(info)
    if not returncol:
        return info
    else:
        return info[returncol]


def download(uuid, slcdir):
    scihub_user, scihub_pass, scihub_url = get_scihub_creds()
    scihub = SentinelAPI(scihub_user, scihub_pass, scihub_url)
    rc = scihub.download(uuid, slcdir)
    return rc


def download_asf(filename, slcdir = '/gws/nopw/j04/nceo_geohazards_vol2/LiCS/temp/SLC', ingest = False):
    # slcdir = os.environ['LiCSAR_SLC']
    wgetpath = os.environ['LiCSARpath']+'/bin/scripts/wget_alaska'
    cmd = 'cd {0}; {1} {2}'.format(slcdir, wgetpath, filename)
    rc = os.system(cmd)
    if ingest:
        filepath = os.path.join(slcdir,filename)
        os.system('arch2DB.py -f '+filepath)
    return rc

def search_alaska(frame, footprint, startdate, enddate, sensType = 'IW'):
    print('performing data discovery using ASF server')
    track = int(frame[0:3])
    trackpre=abs(track-1)
    trackpost=track+1
    strtrack = str(trackpre)+'-'+str(trackpost)
    ascdesc = frame[3]
    if ascdesc == 'D': ascdesc = 'DESCENDING'
    if ascdesc == 'A': ascdesc = 'ASCENDING'
    url = 'https://api.daac.asf.alaska.edu/services/search/param?platform=S1&processingLevel=SLC&output=JSON'
    url = url+'&relativeOrbit='+strtrack
    url = url+'&beamMode='+sensType
    url = url+'&start={0}&end={1}'.format(startdate.strftime('%Y-%m-%d'),enddate.strftime('%Y-%m-%d'))
    url = url + '&flightDirection='+ascdesc
    url = url + '&intersectsWith='+footprint
    r = requests.get(url)
    df = pd.DataFrame.from_dict(r.json()[0])
    return df


def get_images_for_frame(frameName, startdate = dt.datetime.strptime('20141001','%Y%m%d').date(),
             enddate = dt.date.today(), sensType = 'IW', outAspd = False, asf = True):
    #startdate and enddate should be of type datetime.date (but datetime may also work)
    # problem is that only full days are selected, no search by time
    if str(type(startdate)).split("'")[1] == 'str':
        startdate = dt.datetime.strptime(startdate,'%Y%m%d').date()
    if str(type(enddate)).split("'")[1] == 'str':
        enddate = dt.datetime.strptime(enddate,'%Y%m%d').date()
    #one extra date needed for scihub:
    enddate = enddate + dt.timedelta(days=1)
    if enddate > (dt.date.today() - dt.timedelta(days=1)):
        asf = False
        print('checking the latest data - using only scihub')
    #check/update sensType
    if sensType == 'IW':
        if frameName.split('_')[1]=='SM':
            print('the frame is a stripmap')
            sensType = 'SM'
    #nmax = 100
    import LiCSquery as lq
    footprint = lq.get_wkt_boundaries(frameName)
    ascdesc = frameName[3]
    if ascdesc == 'D': ascdesc = 'DESCENDING'
    if ascdesc == 'A': ascdesc = 'ASCENDING'
    track = int(frameName[0:3])
    #startdate = dt.datetime.strptime(startdate,'%Y%m%d').date()
    result = None
    images = None
    track1=track-1
    track2=track+1
    if track2 == 176:
        track2 = 1
    if track1 == 0:
        track1 = 175
    if outAspd or not asf:
        #print('warning, we use outAspd only for EIDP')
        #print('to avoid complications, the search will be performed only through scihub')
        try:
            #use also +-1 track

            scihub_user, scihub_pass, scihub_url = get_scihub_creds()
            scihub = SentinelAPI(scihub_user, scihub_pass, scihub_url)
            result = scihub.query(footprint, date = (startdate.strftime('%Y%m%d'), enddate.strftime('%Y%m%d')), \
                         platformname = 'Sentinel-1', producttype = 'SLC', \
                         relativeorbitnumber = {track1, track, track2}, sensoroperationalmode = sensType, orbitdirection = ascdesc)
            if outAspd:
                images = scihub.to_dataframe(result)
            else:
                df = scihub.to_dataframe(result)
                images = df['title'].values.tolist()
        except:
            print('error in scihub search')
    else:
        try:
            df = search_alaska(frameName, footprint, startdate, enddate, sensType)
            images = df['granuleName'].values.tolist()
        except:
            print('error searching through ASF, trying scihub')
            try:
                scihub_user, scihub_pass, scihub_url = get_scihub_creds()
                scihub = SentinelAPI(scihub_user, scihub_pass, scihub_url)
                result = scihub.query(footprint, date = (startdate.strftime('%Y%m%d'), enddate.strftime('%Y%m%d')), \
                             platformname = 'Sentinel-1', producttype = 'SLC', \
                             relativeorbitnumber = {track1, track, track2}, sensoroperationalmode = sensType, orbitdirection = ascdesc)
                df = scihub.to_dataframe(result)
                images = df['title'].values.tolist()
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
    return images


def get_scihub_creds():
    scihub_url = "https://scihub.copernicus.eu/dhus"
    credentials = os.path.join(os.environ["LiCSAR_configpath"],'scihub_credentials_wget')
    f = open(credentials,'r')
    scihub_user = re.sub('\W+','',f.readline().split('=')[1])
    scihub_pass = re.sub('\W+','',f.readline().split('=')[1])
    f.close()
    return scihub_user, scihub_pass, scihub_url


def get_neodc_path_images(images, file_or_meta = False):
    '''
    work for both one image or image_list
    '''
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
        if file_or_meta:
            if not os.path.exists(neodcpath):
                neodcpath = neodcpath.replace('.zip','.metadata_only.zip')
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
                localfile = os.path.join(os.environ['LiCSAR_SLC'], os.path.basename(imagepath))
                if os.path.exists(localfile):
                    os.system('arch2DB.py -f {}'.format(localfile))
                else:
                    print('Image '+os.path.basename(imagepath)+' is not in neodc, neither downloaded. Including to the list for download')
                    todown.append(os.path.basename(imagepath))
    return todown


def check_and_import_to_licsinfo(frameName, startdate = dt.datetime.strptime('20141001','%Y%m%d').date(), enddate = dt.date.today(), meta = True, reingest = False):
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
    if reingest:
        import framecare as fc
        for fileid in images:
            fc.reingest_file(fileid)
    todown = import_to_licsinfo(images, meta)
    if len(todown)>0:
        print('Summary: There are '+str(len(todown))+' images that were physically acquired but do not exist in CEDA Archive')
    else:
        print('All needed files should exist in CEDA Archive')
    return todown
