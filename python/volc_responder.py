#!/usr/bin/env python

from urllib import response

from lxml import etree
import os
from subprocess import call
import requests
import LiCSquery as lq
from datetime import datetime, timedelta
from batchEnvLib import get_rslc_list
import glob
import pandas as pd

public_path = os.environ['LiCSAR_public']
framelistfile = '/gws/ssde/j25a/nceo_geohazards/vol1/projects/LiCS/volc-proc/active_frames.txt'
customframesfile = '/gws/ssde/j25a/nceo_geohazards/vol1/projects/LiCS/volc-proc/active_frames_custom.txt'
kmzfile = '/gws/ssde/j25a/nceo_geohazards/vol1/projects/LiCS/volc-proc/WeeklyVolcanoGE-Reports.kmz'


import xml.etree.ElementTree as ET

import geopandas as gpd
import pandas as pd
from shapely.geometry import Point

def georss_to_kml(rss_url, output_path="weekly_volcanoes.kml"):
    response = requests.get(rss_url, verify=False, timeout=60)
    response.raise_for_status()
    root = ET.fromstring(response.content)
    ns = {"georss": "http://www.georss.org/georss"}
    rows = []
    for item in root.findall(".//item"):
        point = item.findtext("georss:point", namespaces=ns)
        if not point:
            continue
        lat, lon = map(float, point.split())
        rows.append({
            "name": item.findtext("title", ""),
            "description": item.findtext("description", ""),
            "link": item.findtext("link", ""),
            "geometry": Point(lon, lat),
        })
    gdf = gpd.GeoDataFrame(rows, geometry="geometry", crs="EPSG:4326")
    gdf.to_file(output_path, driver="KML")
    return gdf


def extract_volcanoes_from_RSS(RSS_URL = "https://a.cf11.si.edu/news/WeeklyVolcanoRSS.xml"):
    ''' found by Pedro - thank you'''
    # found in google to use smithosnian data base
    import requests
    import urllib3
    import xml.etree.ElementTree as ET
    import pandas as pd
    import re
    KEYWORDS = [
        "deformation",
        "ground deformation",
        "inflation",
        "deflation",
        "uplift",
        "subsidence",
        "gps",
        "gnss",
        "tilt",
        "geodetic",
        "seismic swarm",
        "swarm",
        "seismicity",
        "seismic",
        "unrest",
        "alert level"
    ]
    urllib3.disable_warnings()
    print("Downloading RSS...")
    response = requests.get(
        RSS_URL,
        verify=False,
        timeout=60
    )
    response.raise_for_status()
    print("Status:", response.status_code)
    root = ET.fromstring(response.content)
    items = list(root.iter("item"))
    print(f"Items found: {len(items)}")
    rows = []
    for item in items:
        title = item.findtext("title", "")
        description = item.findtext("description", "")
        link = item.findtext("link", "")
        volcano = ""
        country = ""
        report_date = ""
        activity_type = ""
        m = re.match(
            r"(.+?)\s+\((.+?)\)\s+-\s+Report for\s+(.*?)\s+-\s+(.*)",
            title
        )
        if m:
            volcano = m.group(1).strip()
            country = m.group(2).strip()
            report_date = m.group(3).strip()
            activity_type = m.group(4).strip()
        text = f"{title} {description}".lower()
        found_keywords = []
        for kw in KEYWORDS:
            if kw.lower() in text:
                found_keywords.append(kw)
        rows.append({
            "volcano": volcano,
            "country": country,
            "report_date": report_date,
            "activity_type": activity_type,
            "keywords": ";".join(found_keywords),
            "n_keywords": len(found_keywords),
            "link": link
        })
    df = pd.DataFrame(rows)
    df.sort_values(
        by="n_keywords",
        ascending=False,
        inplace=True
    )
    return df


def download_volc_kml(kmzfile = 'WeeklyVolcanoGE-Reports.kmz'):
    #thanks Fabien for this!
    url = 'https://volcano.si.edu/news/'+kmzfile
    kmlfile = kmzfile.replace('.kmz','.kml')
    if os.path.exists(kmzfile):
        os.remove(kmzfile)
    if os.path.exists(kmlfile):
        os.remove(kmlfile)
    #if use_playwright:
    #    from playwright.sync_api import sync_playwright
    #    with sync_playwright() as p:
    # bot protection... need more complex:
    headers = {
        "User-Agent": (
            "Mozilla/5.0 (Windows NT 10.0; Win64; x64) "
            "AppleWebKit/537.36 (KHTML, like Gecko) "
            "Chrome/139.0.0.0 Safari/537.36"
        ),
        "Referer": "https://volcano.si.edu/",
    }
    r = requests.get(url, headers=headers, allow_redirects=True)
    #print(r.status_code)
    #print(r.headers.get("content-type"))
    with open(kmzfile, "wb") as f:
        f.write(r.content)
    # request = requests.get(url)
    if r.status_code == 200:
        # os.system('wget {0} >/dev/null 2>/dev/null'.format(url))
        os.system('unzip {0} >/dev/null 2>/dev/null'.format(kmzfile))
        os.system('rm {0}'.format(kmzfile))
        os.system('iconv -f iso-8859-1 -t utf-8 {0} > {1}'.format(kmlfile, kmlfile+'.utf8'))
        os.system('rm {0}'.format(kmlfile))
        os.rename(kmlfile+'.utf8', kmlfile)
    else:
        print('some URL problem - not downloading')
        print('error: '+str(r))
        return False
    return True


def get_frames_from_kml(filename = 'WeeklyVolcanoGE-Reports.kml'):
    doc = etree.parse(filename)
    nmsp = '{http://www.opengis.net/kml/2.2}'
    frames_set = set()
    #
    for pm in doc.iterfind('.//{0}Placemark'.format(nmsp)):
        for ls in pm.iterfind('{0}Point/{0}coordinates'.format(nmsp)):
            coordinate = ls.text.strip().replace('\n','')
            lon = float(coordinate.split(',')[0])
            lat = float(coordinate.split(',')[1])
            #get frames for this lon lat:
            frames = lq.get_frames_in_polygon(lon, lon+0.1, lat, lat+0.1)
            frames = lq.sqlout2list(frames)
            for frame in frames:
                frames_set.add(frame)
    return list(frames_set)


def update_framelist_fabien(volcdir = '/gws/ssde/j25a/nceo_geohazards/vol1/projects/LiCS/volc-proc/list_database'):
    '''
    this will update list of frames that will be used for given volcano in the volcano portal
    '''
    volcfiles = glob.glob(os.path.join(volcdir,'*volcano*txt'))
    framefiles = []
    for volc in volcfiles:
        framef = volc.replace('volcano','frame')
        print('updating '+framef)
        framefiles.append(framef)
        if os.path.exists(framef):
            os.remove(framef)
        volcs = pd.read_csv(volc,sep=r'\s+',header=None)
        allframes = []
        for i,v in volcs.iterrows():
            lat=v[1]
            lon=v[2]
            frames = lq.get_frames_in_lonlat(lon,lat)
            frames = lq.sqlout2list(frames)
            allframes = allframes + frames
        allframes=list(set(allframes))
        for frame in allframes:
            rc = os.system('echo {0} >> {1}'.format(frame,framef))


def get_indate(frame, numepochs = 3):
    #this is to get indate for licsar_make_frame - either 90 days ago
    #or day before the last three existing RSLCs
    try:
        rslcs = get_rslc_list(frame, True)
    except:
        rslcs = ''
    #checking if TOGETHER WITH MASTER there is more than numepochs epochs
    if len(rslcs) > numepochs:
        indate = rslcs[-numepochs]
        indate = datetime.strptime(indate, '%Y%m%d') - timedelta(days=1)
    else:
        indate = datetime.now()-timedelta(days=90)
    return indate.date()

def main():
    print('getting Weekly Active Volcanoes from Smithonian RSS feed')
    try:
        kmlfile = 'WeeklyVolcanoGE-fromRSS.kml'
        gdf = georss_to_kml(rss_url = "https://a.cf11.si.edu/news/WeeklyVolcanoRSS.xml", output_path=kmlfile)
        frames = get_frames_from_kml(kmlfile)
    except:
        print('some error getting the RSS feed - trying to download the kmz file instead')
        if not download_volc_kml(kmzfile):
            print('some error in downloading weekly kmz - using the previous one')    
        print('getting list of frames covering the active volcanoes')
        kmlfile = kmzfile.replace('.kmz','.kml')
        frames = get_frames_from_kml(kmlfile)
    print('There are {0} frames to process over the active volcanoes'.format(str(len(frames))))
    with open(customframesfile, 'r') as f:
        custom_frames = f.read().splitlines()
    for cf in custom_frames:
        if cf in frames:
            custom_frames.remove(cf)
    print('adding other {} custom frames to the set'.format(len(custom_frames)))
    frames = frames + custom_frames
    print('The final list is saved to: {}'.format(framelistfile))
    with open(framelistfile, 'w') as f:
        for frame in frames:
            f.write("%s\n" % frame)
    #set processing dates - last day is 'tomorrow', first day depends on last two images of frames
    offdate = datetime.now()+timedelta(days=1)
    for frame in frames:
        track = str(int(frame[0:3]))
        if not os.path.exists(os.path.join(public_path,track,frame)):
            print('Frame '+frame+' was (probably) not initiated. please try init it first using:') #, trying to do it automatically')
            cmd = 'licsar_initiate_new_frame.sh {0}'.format(frame)
            print(cmd)
            #os.system(cmd)
        if os.path.exists(os.path.join(public_path,track,frame)):
            indate = get_indate(frame)
            print('..preparing frame {0} and sending processing jobs to LOTUS'.format(frame))
            # we used -P before, now -R for comet_ ...responder
            #try:
            #    volcpath = os.path.join(os.environ['LiCSAR_procdir'],'..','..','projects','LiCS','volc-proc')
            #except:
            #    volcpath = '.'
            volcpath = os.path.join(os.environ['LiCSAR_procdir'],'..','..','volc-proc')
            cmd = 'licsar_make_frame.sh -S -N -f {0} 0 1 {1} {2} > {3}/volcresp.{0}.log 2> {3}/volcresp.{0}.err'.format(frame,str(indate),str(offdate.date()), volcpath)
            print(cmd)
            os.system(cmd)
    #first update the frame lists - only once per month...:
    nowis = datetime.now()
    if nowis.day == 9 or nowis.day == 29:
        print('updating frame list first')
        update_framelist_fabien()
    # but update also tifs for volc portal, using Fabien's script, every 10 days
    #if nowis.day in [1,8,15,22,29]:
    #if nowis.day in [5,15,25]:
    #if nowis.day in [5,15,24, 26, 28]:
    #print('done. now updating the png files - Fabien script')
    #cmd = 'cd /gws/ssde/j25a/nceo_geohazards/vol1/projects/LiCS/volc-proc/python_script/FINAL_scripts; sbatch run_LiCSAR_volcano_makefigure_final_withGACOS.LOTUS2.sh'
    #os.system(cmd)
    print('ok, everything finished')


if __name__ == '__main__':
    main()
