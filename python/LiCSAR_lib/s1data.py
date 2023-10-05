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



def get_bperps_asf(product_id):
    try:
        import asf_search as asf
        mdate=product_id.split('T')[0].split('_')[-1]
        print('searching for stack through ASF')
        reference = asf.product_search(product_id+'-SLC')[0]
        stack = reference.stack()
        bls=asf.get_baseline_from_stack(reference, stack)
        bperps=[]
        btemps=[]
        bdates=[]
        mdates=[]
        for prod in bls[0]:
            bperp=prod.properties.get('perpendicularBaseline')
            btemp=prod.properties.get('temporalBaseline')
            bdate=prod.properties.get('fileID').split('T')[0].split('_')[-1]
            btemps.append(btemp)
            bperps.append(bperp)
            bdates.append(bdate)
            mdates.append(mdate)
        pdict = {'ref_date': mdates, 'date': bdates, 'bperp': bperps, 'btemp': btemps}
        bperpd = pd.DataFrame(pdict)
        return bperpd
    except:
        print('ERROR: perhaps install asf_search, or it did not find the data')
        return False
    


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


def get_epochs_for_frame(frame, startdate = dt.datetime.strptime('20141001','%Y%m%d').date(), enddate = dt.date.today(), returnAsDate = False):
    new_images = get_images_for_frame(frame, startdate, enddate)
    epochs = [i.split('_')[5].split('T')[0] for i in new_images]
    epochs = list(set(epochs))
    if returnAsDate:
        return [dt.date(int(a[:4]),int(a[4:6]),int(a[6:8])) for a in epochs]
    else:
        return epochs


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


# addition to CDSE - based on https://github.com/DHI-GRAS/creodias-finder/blob/main/creodias_finder
import requests
from tqdm import tqdm
import shutil
from pathlib import Path
import concurrent.futures
import datetime
from six.moves.urllib.parse import urlencode
from six import string_types
import dateutil.parser
from shapely.geometry import shape
import re

API_URL = (
    "https://catalogue.dataspace.copernicus.eu/resto/api/collections/{collection}"
    "/search.json?maxRecords=1000"
)
ONLINE_STATUS_CODES = "34|37|0"
DOWNLOAD_URL = "https://zipper.creodias.eu/download"
TOKEN_URL = "https://identity.dataspace.copernicus.eu/auth/realms/CDSE/protocol/openid-connect/token"


def cdse_query(
    collection,
    start_date=None,
    end_date=None,
    geometry=None,
    status=ONLINE_STATUS_CODES,
    **kwargs ):
    """Query the EOData Finder API
    
    Parameters
    ----------
    collection: str, optional
        the data collection, corresponding to various satellites
    start_date: str or datetime
        the start date of the observations, either in iso formatted string or datetime object
    end_date: str or datetime
        the end date of the observations, either in iso formatted string or datetime object
        if no time is specified, time 23:59:59 is added.
    geometry: WKT polygon or object impementing __geo_interface__
        area of interest as well-known text string
    status : str
        allowed online/offline statuses (|-separated for OR)
    **kwargs
        Additional arguments can be used to specify other query parameters,
        e.g. productType=L1GT
        See https://creodias.eu/eo-data-finder-api-manual for a full list
    
    Returns
    -------
    dict[string, dict]
        Products returned by the query as a dictionary with the product ID as the key and
        the product's attributes (a dictionary) as the value.
    """
    query_url = _build_query(
        API_URL.format(collection=collection),
        start_date,
        end_date,
        geometry,
        status,
        **kwargs,
    )
    
    query_response = {}
    while query_url:
        response = requests.get(query_url)
        response.raise_for_status()
        data = response.json()
        for feature in data["features"]:
            query_response[feature["id"]] = feature
        query_url = _get_next_page(data["properties"]["links"])
    return query_response


def _build_query(
    base_url, start_date=None, end_date=None, geometry=None, status=None, **kwargs):
    query_params = {}
    
    if start_date is not None:
        start_date = _parse_date(start_date)
        query_params["startDate"] = start_date.isoformat()
    if end_date is not None:
        end_date = _parse_date(end_date)
        end_date = _add_time(end_date)
        query_params["completionDate"] = end_date.isoformat()
    
    if geometry is not None:
        query_params["geometry"] = _parse_geometry(geometry)
    
    if status is not None:
        query_params["status"] = status
    
    for attr, value in sorted(kwargs.items()):
        value = _parse_argvalue(value)
        query_params[attr] = value
    
    url = base_url
    if query_params:
        url += f"&{urlencode(query_params)}"
    
    return url


def _get_next_page(links):
    for link in links:
        if link["rel"] == "next":
            return link["href"]
    return False


def _parse_date(date):
    if isinstance(date, datetime.datetime):
        return date
    elif isinstance(date, datetime.date):
        return datetime.datetime.combine(date, datetime.time())
    try:
        return dateutil.parser.parse(date)
    except ValueError:
        raise ValueError(
            "Date {date} is not in a valid format. Use Datetime object or iso string"
        )


def _add_time(date):
    if date.hour == 0 and date.minute == 0 and date.second == 0:
        date = date + datetime.timedelta(hours=23, minutes=59, seconds=59)
        return date
    return date


def _tastes_like_wkt_polygon(geometry):
    try:
        return geometry.replace(", ", ",").replace(" ", "", 1).replace(" ", "+")
    except Exception:
        raise ValueError("Geometry must be in well-known text format")


def _parse_geometry(geom):
    try:
        # If geom has a __geo_interface__
        return shape(geom).wkt
    except AttributeError:
        if _tastes_like_wkt_polygon(geom):
            return geom
        raise ValueError(
            "geometry must be a WKT polygon str or have a __geo_interface__"
        )


def _parse_argvalue(value):
    if isinstance(value, string_types):
        value = value.strip()
        if not any(
            value.startswith(s[0]) and value.endswith(s[1])
            for s in ["[]", "{}", "//", "()"]
        ):
            value.replace(" ", "+")
        return value
    elif isinstance(value, (list, tuple)):
        # Handle value ranges
        if len(value) == 2:
            value = "[{},{}]".format(*value)
            return value
        else:
            raise ValueError(
                "Invalid number of elements in list. Expected 2, received "
                "{}".format(len(value))
            )
    else:
        raise ValueError(
            "Additional arguments can be either string or tuple/list of 2 values"
        )

def _get_token(username, password):
    token_data = {
        "client_id": "cdse-public",
        "username": username,
        "password": password,
        "grant_type": "password",
    }
    response = requests.post(TOKEN_URL, data=token_data).json()
    try:
        return response["access_token"]
    except KeyError:
        raise RuntimeError(f"Unable to get token. Response was {response}")


def cdse_download(uid, username, password, outfile, show_progress=True, token=None):
    """Download a file from CreoDIAS to the given location

    Parameters
    ----------
    uid:
        CreoDIAS UID to download
    username:
        Username
    password:
        Password
    outfile:
        Path where incomplete downloads are stored
    """
    token = token if token else _get_token(username, password)
    url = f"{DOWNLOAD_URL}/{uid}?token={token}"
    _download_raw_data(url, outfile, show_progress)


def cdse_download_list(uids, username, password, outdir, threads=1, show_progress=True):
    """Downloads a list of UIDS
    
    Parameters
    ----------
    uids:
        A list of UIDs
    username:
        Username
    password:
        Password
    outdir:
        Output direcotry
    threads:
        Number of simultaneous downloads
    
    Returns
    -------
    dict
        mapping uids to paths to downloaded files
    """
    if show_progress:
        pbar = tqdm(total=len(uids), unit="files")
    
    token = _get_token(username, password)
    
    def _download(uid):
        outfile = Path(outdir) / f"{uid}.zip"
        download(
            uid, username, password, outfile=outfile, show_progress=False, token=token
        )
        if show_progress:
            pbar.update(1)
        return uid, outfile
    
    with concurrent.futures.ThreadPoolExecutor(threads) as executor:
        paths = dict(executor.map(_download, uids))
    
    return paths


def _download_raw_data(url, outfile, show_progress):
    """Downloads data from url to outfile.incomplete and then moves to outfile"""
    outfile_temp = str(outfile) + ".incomplete"
    try:
        downloaded_bytes = 0
        with requests.get(url, stream=True, timeout=100) as req:
            print(req.status_code)
            with tqdm(unit="B", unit_scale=True, disable=not show_progress) as progress:
                chunk_size = 2**20  # download in 1 MB chunks
                with open(outfile_temp, "wb") as fout:
                    for chunk in req.iter_content(chunk_size=chunk_size):
                        if chunk:  # filter out keep-alive new chunks
                            fout.write(chunk)
                            progress.update(len(chunk))
                            downloaded_bytes += len(chunk)
        shutil.move(outfile_temp, outfile)
    finally:
        try:
            Path(outfile_temp).unlink()
        except OSError:
            pass


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
