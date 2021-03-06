#!/usr/bin/env python

import numpy as np
import os, re
import argparse
import datetime
from html.parser import *

import requests
# Add this line because insecure requests warning
# suggested solution in http://stackoverflow.com/questions/27981545/surpress-insecurerequestwarning-unverified-https-request-is-being-made-in-pytho
requests.packages.urllib3.disable_warnings()

# Old files seem not be found with this script

maxpages = 10 #will search through a maximum of 10 pages (if more expected, up this number)
server = 'https://qc.sentinel1.copernicus.eu/'
orbitMap = [('precise','aux_poeorb'),('restituted','aux_resorb')]
datefmt = "%Y%m%dT%H%M%S"
queryfmt = "%Y-%m-%d"

def cmdLineParse():
    '''
    Command line parser.
    '''
    parser = argparse.ArgumentParser(description='Fetch orbits corresponding to given SAFE package')
    parser.add_argument('-i', '--input', dest='input', type=str, required=True,
            help='Path to SAFE package of interest')
    parser.add_argument('-o', '--output', dest='outdir', type=str, default='.',
            help='Path to output directory')
    return parser.parse_args()

def FileToTimeStamp(safename):
    '''
    Return timestamp from SAFE name.
    '''
    fields = safename.split('_')
    tstamp = datetime.datetime.strptime(fields[-4], datefmt)
    return tstamp

class MyHTMLParser(HTMLParser):
    def __init__(self):
        HTMLParser.__init__(self)
        self.fileList = []
        self.in_td = False
        self.in_a = False

    def handle_starttag(self, tag, attrs):
        if tag == 'td':
            self.in_td = True
        elif tag == 'a':
            self.in_a = True

    def handle_data(self,data):
        if self.in_td and self.in_a:
            if re.search(r"S1._OPER",data):
                self.fileList.append(data.strip())

    def handle_tag(self, tag):
        if tag == 'td':
            self.in_td = False
            self.in_a = False
        elif tag == 'a':
            self.in_a = False

def download_file(url, outdir='.', session=None):
    '''
    Download file to specified directory.
    '''
    if session is None:
        session = requests.session()

    path = os.path.join(outdir, os.path.basename(url))
    print(('Downloading URL: ', url))
    request = session.get(url, stream=True, verify=False)

    try:
        val = request.raise_for_status()
        success = True
    except:
        success = False

    if success:
        with open(path,'wb') as f:
            for chunk in request.iter_content(chunk_size=1024):
                if chunk:
                    f.write(chunk)
                    f.flush()

    return success

if __name__ == '__main__':
    '''
    Main driver.
    '''
    inps = cmdLineParse()
    sat = inps.input[:3]
    fileTS = FileToTimeStamp(inps.input)
    print(('Reference time: ', fileTS))
    match = None
    session = requests.Session()
    for spec in orbitMap:
        oType = spec[0]
        if oType == 'precise':
            delta = datetime.timedelta(days=2)
        elif oType == 'restituted':
            delta = datetime.timedelta(hours=3)

        timebef = (fileTS - delta).strftime(queryfmt)
        timeaft = (fileTS + delta).strftime(queryfmt)
        url = server + spec[1]
        
        for page in range (1,maxpages):
            query = url + '/?page=' + str(page) + '&mission='+sat+'&validity_start_time={0}..{1}'.format(timebef, timeaft)
            print(query)
            success = False
            match = None
            try:
                print(('Querying for {0} orbits'.format(oType)))
                r = session.get(query, verify=False)
                r.raise_for_status()
                res = r.text
                parser = MyHTMLParser()
                parser.feed(r.text)
                match = None
                for result in parser.fileList:
                    print(result)
                    fields = result.split('_')
                    taft = datetime.datetime.strptime(fields[-1][0:15], datefmt)
                    tbef = datetime.datetime.strptime(fields[-2][1:16], datefmt)
                    if (tbef <= fileTS) and (taft >= fileTS):
                        match = os.path.join(url, result)
                        break

                if match is not None:
                    success = True
                    break

            except :
                pass

        if match is None:
            print(('Failed to find {0} orbits for Time {1}'.format(oType, fileTS)))

        if success:
            break

    if match:
        res = download_file(match, inps.outdir, session=session)

        if res is False:
            print(('Failed to download URL: ', match))

    session.close()
