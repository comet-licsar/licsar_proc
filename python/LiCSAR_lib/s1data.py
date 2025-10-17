#!/usr/bin/env python

import os
import datetime as dt
# from sentinelsat.sentinel import SentinelAPI # this was beautiful toolbox, sadly changed to CDSE..
import re
import pandas as pd
import requests
#import arch2DB
#here is the nostdout function:
from LiCSAR_misc import *
import asf_search as asf

# need to update this one
def get_info_pd(fileid = 'S1A_IW_SLC__1SDV_20210908T235238_20210908T235305_039597_04AE3C_4CA7', returncol = None):
    filename = fileid+'.SAFE'
    #scihub_user, scihub_pass, scihub_url = get_scihub_creds()
    #scihub = SentinelAPI(scihub_user, scihub_pass, scihub_url)
    #info = scihub.query(filename = filename)
    #info = scihub.to_dataframe(info)
    json = requests.get("https://catalogue.dataspace.copernicus.eu/odata/v1/Products?$filter=Name eq '"+filename+"'").json()
    info = pd.DataFrame.from_dict(json["value"]).head(1)
    # CDSE e.g. "geography'SRID=4326;MULTIPOLYGON (((93.754837 35.530998, 94.143021 37.151909, 91.328789 37.546482, 91.001236 35.9272, 93.754837 35.530998)))'"
    info['footprint']=info.Footprint.values[0].split(';')[-1][:-1]
    if not returncol:
        return info
    else:
        return info[returncol]

'''
# SciHub is gone - sentinelsat not updated
def download(uuid, slcdir):
    scihub_user, scihub_pass, scihub_url = get_scihub_creds()
    scihub = SentinelAPI(scihub_user, scihub_pass, scihub_url)
    rc = scihub.download(uuid, slcdir)
    return rc
'''

def download(filename, slcdir = '/gws/nopw/j04/nceo_geohazards_vol2/LiCS/temp/SLC', ingest = False, provider='cdse'):
    '''wrapper to wget commands. the provider must be one of ['cdse', 'alaska']
    '''
    # slcdir = os.environ['LiCSAR_SLC']
    wgetpath = os.environ['LiCSARpath']+'/bin/scripts/wget_'+provider
    cmd = 'cd {0}; {1} {2}'.format(slcdir, wgetpath, filename)
    rc = os.system(cmd)
    filepath = os.path.join(slcdir,filename)
    if ingest:
        os.system('arch2DB.py -f '+filepath)
    return filepath


def get_bperps_asf(product_id):
    try:
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
        bperpd=bperpd.dropna()
        bperpd=bperpd.drop_duplicates()
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
    tracks = [trackpre, track, trackpost]
    #strtrack = str(trackpre)+'-'+str(trackpost)
    ascdesc = frame[3]
    if ascdesc == 'D': ascdesc = 'DESCENDING'
    if ascdesc == 'A': ascdesc = 'ASCENDING'
    #url = 'https://api.daac.asf.alaska.edu/services/search/param?platform=S1&processingLevel=SLC&output=JSON'
    #url = url+'&relativeOrbit='+strtrack
    #url = url+'&beamMode='+sensType
    #url = url+'&start={0}&end={1}'.format(startdate.strftime('%Y-%m-%d'),enddate.strftime('%Y-%m-%d'))
    #url = url + '&flightDirection='+ascdesc
    #url = url + '&intersectsWith='+footprint
    #r = requests.get(url)
    # or using asf search:
    r = asf.geo_search(relativeOrbit=tracks, beamMode=sensType,
                   start = startdate, end = enddate,
                   flightDirection=ascdesc, intersectsWith=footprint,
                   platform='S1', processingLevel='SLC')
    images = []
    for gg in r:
        images.append(gg.properties['sceneName'])
    #df = pd.DataFrame.from_dict(r.json()[0])
    #return df
    return images


def get_epochs_for_frame(frame, startdate = dt.datetime.strptime('20141001','%Y%m%d').date(), enddate = dt.date.today(), returnAsDate = False):
    new_images = get_images_for_frame(frame, startdate, enddate)
    epochs = [i.split('_')[5].split('T')[0] for i in new_images]
    epochs = list(set(epochs))
    if returnAsDate:
        return [dt.date(int(a[:4]),int(a[4:6]),int(a[6:8])) for a in epochs]
    else:
        return epochs


def get_images_for_footprint(frameName, footprint, startdate = dt.datetime.strptime('20141001','%Y%m%d').date(),
                         enddate = dt.date.today(), sensType = 'IW'):
    '''frameName can be a fake one, e.g. '018D' is enough. footprint is a POLYGON WKT, e.g.
    bidsgpd = fc.bursts2geopandas([burstid])
    footprint = bidsgpd.geometry[0].wkt
    '''
    images = search_alaska(frameName, footprint, startdate, enddate, sensType)
    return images


def get_images_for_burst(bidtanx, orbdir = 'A'):
    ''' orbdir for either A or D (scending) orbit direction'''
    import framecare as fc
    bidsgpd = fc.bursts2geopandas([bidtanx])
    footprint = bidsgpd.geometry[0].wkt
    relorb = bidtanx.split('_')[0]
    if int(relorb)<100:
        relorb = '0'+relorb
    if int(relorb) < 10:
        relorb = '0' + relorb
    relorb = relorb+orbdir
    return get_images_for_footprint(relorb, footprint)


def get_images_for_frame(frameName, startdate = dt.datetime.strptime('20141001','%Y%m%d').date(),
             enddate = dt.date.today(), sensType = 'IW', outAspd = False, asf = True):
    ''' Will get filenames from CDSE and ASF for the frame. Note this is based on frame polygon, overlapping bursts might cause issues!

    Args:
        frameName: str
        startdate: dt.date
        enddate: dt.date
        sensType: 'IW' or 'SM'
        outAspd: if True, it will use only CDSE (skip ASF) and return the outputs as more complete pd.DataFrame
        asf: would use ONLY ASF (will avoid CDSE) for the search (does not include everything...), otherwise CDSE+ASF (or purely CDSE in case of outAspd==True)

    Returns:
        list or pd.DataFrame
    '''
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
        print('checking the latest data - using only CDSE')
        #print('SHOULD DO FROM CDSE - keeping ASF for now, no data for last 24(or 48?) hours') 
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
        # 2023-11-06 - searching through CDSE, solution by Manu Delgado Blasco (many thanks!) https://github.com/sentinelsat/sentinelsat/issues/583
        # cannot find docs for OData! e.g. - use of 'sensoroperationalMode' ends by Invalid field: sensoroperationalMode (BUT WHAT ARE VALID FIELDS????)
        topp = 400 # max 1000   ### 2023-11-16 fix
        cdsequery = f"https://catalogue.dataspace.copernicus.eu/odata/v1/Products?$filter=Collection/Name eq 'SENTINEL-1' and " \
        "OData.CSC.Intersects(area=geography'SRID=4326;{0}') and ContentDate/Start gt {1}T00:00:00.000Z and ContentDate/Start lt {2}T00:00:00.000Z " \
        "and Attributes/OData.CSC.StringAttribute/any(att:att/Name eq 'productType' and att/OData.CSC.StringAttribute/Value eq 'SLC') " \
        "and (Attributes/OData.CSC.StringAttribute/any(att:att/Name eq 'relativeOrbitNumber' and att/OData.CSC.IntegerAttribute/Value eq {3})" \
        " or Attributes/OData.CSC.StringAttribute/any(att:att/Name eq 'relativeOrbitNumber' and att/OData.CSC.IntegerAttribute/Value eq {4})" \
        " or Attributes/OData.CSC.StringAttribute/any(att:att/Name eq 'relativeOrbitNumber' and att/OData.CSC.IntegerAttribute/Value eq {5})) " \
        "and Attributes/OData.CSC.StringAttribute/any(att:att/Name eq 'orbitDirection' and att/OData.CSC.StringAttribute/Value eq '{6}') " \
        "and Attributes/OData.CSC.StringAttribute/any(att:att/Name eq 'operationalMode' and att/OData.CSC.StringAttribute/Value eq '{7}')" \
        "&$top={8}".format(footprint, str(startdate), str(enddate),
            str(track1),str(track),str(track2),
            ascdesc, sensType,
            str(topp)
            )
        json = requests.get(cdsequery).json()
        # for list of params in the new CDSE (THEY ARE NOT PUBLISHED!!!!!!!!!!!!!!!!! in NOV 2023 when SciHub is deactivated!!!!! this is not nice approach from ESA)
        # can be found after 'expanding the metadata', e.g. here:
        # https://catalogue.dataspace.copernicus.eu/odata/v1/Products?%24filter=contains(Name,%27S1A_EW_GRD%27)%20and%20ContentDate/Start%20gt%202022-05-03T00:00:00.000Z%20and%20ContentDate/Start%20lt%202022-05-03T12:00:00.000Z&%24expand=Attributes
        #
        dframe = pd.DataFrame.from_dict(json["value"])
        if dframe.empty:
            print('CDSE: empty output')
            return False
        i = 0
        dframefull = dframe.copy()
        while not dframe.empty:
            i = i+1
            json = requests.get(cdsequery+"&$skip="+str(i*topp)).json()
            dframe = pd.DataFrame.from_dict(json["value"])
            dframefull = pd.concat([dframefull, dframe], ignore_index=True)
        #
        dframefull['title'] = dframefull['Name'].apply(lambda x: x.split('.')[0])
        if outAspd:
            return dframefull
        else:
            images = dframefull['title'].values.tolist()
            # DEBUG: ASF uses a bit different filename. So adding this here, as ASF is used as backup (and I don't know how to search with filename from ASF to do it through wget_alaska.sh
            try:
                print('CDSE search complete, adding also from ASF (as some filenames there differ in the last 4 digits)')
                images += search_alaska(frameName, footprint, startdate, enddate, sensType)
                #images += df['granuleName'].values.tolist()
            except:
                print('error in connection to ASF')
            return list(set(images))
        '''
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
        '''
    else:
        try:
            images = search_alaska(frameName, footprint, startdate, enddate, sensType)
            #images = df['granuleName'].values.tolist()
        except:
            print('error searching through ASF, cancelling') #', trying scihub')
            '''
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
            '''
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

'''
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
'''

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


def import_to_licsinfo(images, meta = True, extradirs = [os.environ['LiCSAR_SLC']]): #,'/work/xfc/vol5/user_cache/earmla/SLC']):
    ''' if meta, it will import either real zip or metafile (if the zip is not existing)'''
    # updating extradirs
    try:
        extradir=os.path.join(os.environ['XFCPATH'], 'SLC')
        if os.path.exists(extradir):
            extradirs.append(extradir)
        #extradirs2 = []
        #efile = os.path.join(os.environ['LiCSAR_configpath'],'autodownloaddirs')
        #with open(efile) as f:
        #    line = f.readline().split()
        #    extradirs2.append(line[0])
        #extradirs = extradirs+extradirs2
    except:
        print('')
        #print('error reading extra dirs from '+efile)
    extradirs = list(set(extradirs))
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
                for otherdir in extradirs:
                    localfile = os.path.join(otherdir, os.path.basename(imagepath))
                    if os.path.exists(localfile):
                        break
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
    images = get_images_for_frame(frameName, startdate, enddate, sensType, asf=False)  # I need to extract everything from CDSE(+ASF)
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
    todown = import_to_licsinfo(images, meta)  # [ 'S1A...zip', ...]
    if todown:
        import LiCSquery as lq
        print('Ensuring only images with frame bursts are returned')
        framebursts = lq.sqlout2list(lq.get_bidtanxs_in_frame(frameName))
        todowncp = todown.copy()
        for im in todowncp:
            bursts = lq.sqlout2list(lq.get_bursts_in_file(im))
            if not bursts:
                print('error checking file '+im)
            else:
                isinframe = False
                for b in bursts:
                    if b in framebursts:
                        isinframe = True
                        break
                if not isinframe:
                    todown.remove(im)
    if len(todown)>0:
        print('Summary: There are '+str(len(todown))+' images that were physically acquired but do not exist in CEDA Archive')
    else:
        print('All needed files should exist in CEDA Archive')
    return todown
