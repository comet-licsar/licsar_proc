#!/usr/bin/env python

# This file contains the definition of the SLCzip class

import zipfile
import os
import numpy as np
import itertools
import mgrs
import pdb
try:
    import xml.etree.cElementTree as et
except ImportError:
    import xml.etree.ElementTree as et
from math import *

"""
Version 0.1    01 Feb 2016    ELH    Test release, class defined for DB ingestion - much is based on code by Karsten Spaans
Version 0.2    21 Apr 2016    ELH    Added significant number of new functions to original code.
Version 0.3    26 Apr 2016    ELH    Added function to get IPF version
Version 0.4    21 Jun 2016    ELH    Catches errors due to corrupt files
Version 0.5    07 Oct 2016    ELH    Fix for incorrect scene centres over dateline
Version 0.6    08 Nov 2016    ELH    Fix for S1B relative orbit number
"""

class SLCzip(object):
    # This defines the SLC zip file and methods for extracting its data.
    
    def __init__(self, slcfile):
        self.slcfile = slcfile
        #pdb.set_trace()
        try:
            self.manifest = self.extract_manifest()
        except:
            raise
        try:
            self.annotation = self.extract_annotation()
        except:
            raise

    def get_namebase(self):
        # Get the base of the filename, we do this twice as sometimes we have a .SAFE.zip, sometimes, just a .zip
        name = os.path.splitext(self.slcfile)[0]
        namebase = os.path.splitext(name)[0]
        basenm = os.path.basename(namebase)
        return basenm        

    def extract_manifest(self):
        # Extracts the manifest file to memory
        namebase = self.get_namebase()
        manifest = namebase + '.SAFE/manifest.safe'
        try:
            zf = zipfile.ZipFile(self.slcfile)
        except zipfile.BadZipfile as e:
            print('ERROR: Corrupt zipfile: {0}'.format(namebase))
            raise
        try:
            data = zf.read(manifest)
            #print data
            return data
        except KeyError as e:
            print('ERROR: Did not find %s in zipfile' % namebase)
            print('    {}'.format(e.message))
            raise 
        finally:
            zf.close()

    def list_annotation(self):
        # lists the annotation filenames in the archive in an ordered list
        # The order is preserved throughout the Class.
        namebase = self.get_namebase()
        try:
            zf = zipfile.ZipFile(self.slcfile)
        except zipfile.BadZipfile as e:
            print('ERROR: Corrupt zipfile: {0}'.format(namebase))
            raise 
        content = zf.namelist()
        zf.close()
        ann_list_copol = []
        ann_list_crosspol = []
        ann_list = []
        for item in content:
            if 'annotation/s1' in item:
                if ('vv' in item) or ('hh' in item):
                    ann_list_copol.append(item)
                else:
                    ann_list_crosspol.append(item)
        ann_list_copol.sort()
        ann_list_crosspol.sort()
        ann_list = ann_list_copol + ann_list_crosspol
        return ann_list

    def extract_annotation(self):
        # extracts annotation files and returns a list with the XML content of each
        namebase = self.get_namebase()
        ann_list = self.list_annotation()
        try:
            zf = zipfile.ZipFile(self.slcfile)
        except zipfile.BadZipfile as e:
            print('ERROR: Corrupt zipfile: {0}'.format(namebase))
            raise
        data_list = []
        for ann in ann_list:
            try:
                data = zf.read(ann)
                data_list.append(data)
            except KeyError as e:
                print('ERROR: Did not find %s in zipfile' % namebase)
                print('    {}'.format(e.message))
        zf.close()
        return data_list

    def get_orbit(self):
        # returns the relative orbit calculated from the orbit in the manifest file
        root = et.fromstring(self.manifest)
        metadatasection = root.find('metadataSection')
        name = self.get_namebase()
        if 'S1A' in name:
            for el in metadatasection:
                if el.attrib['ID'] == 'measurementOrbitReference':
                    for orbel in el[0][0][0].iter():
                        if 'orbitNumber' in orbel.tag and orbel.attrib['type'] == 'stop':
                            orbnumber = int(orbel.text)
                            cycno = int(orbnumber/175)
                            natorb = (int(orbnumber)-175*cycno)-72
                            if natorb <= 0:
                                relorb = natorb + 175
                            else:
                                relorb = natorb
        elif 'S1B' in name:
            for el in metadatasection:
                if el.attrib['ID'] == 'measurementOrbitReference':
                    for orbel in el[0][0][0].iter():
                        if 'relativeOrbitNumber' in orbel.tag and orbel.attrib['type'] == 'stop':
                            relorb = int(orbel.text)
        else:
            return -1
        return relorb

    def get_polid(self):
        # returns the polarisation of each annotation file
        polid_list = []
        for item in self.annotation:
            root = et.fromstring(item)
            polid_list.append(root.find('adsHeader').find('polarisation').text)
        return polid_list
    
    def get_swathid(self):
        # returns the swath of each annotation file
        swathid_list = []
        for item in self.annotation:
            root = et.fromstring(item)
            swathid_list.append(root.find('adsHeader').find('swath').text)
        return swathid_list
    
    def get_orbdir(self):
        # returns the orbit direction from the manifest file
        root = et.fromstring(self.manifest)
        for elem in root.iter('metadataObject'):
            if elem.attrib['ID'] == 'measurementOrbitReference':
                for orbel in elem[0][0][0].iter():
                    if 'pass' in orbel.tag:
                        orbdir = (orbel.text)[0]
        return orbdir
            
    def get_sensdate(self):
        # returns a list of sensing dates from each annotation file
        sensdate_list = []
        s = (' ')
        for item in self.annotation:
            root = et.fromstring(item)
            sensdate = root.find('adsHeader').find('startTime').text[:10]
            #print sensdate
            senstime = root.find('adsHeader').find('startTime').text[11:]
            sensdate_list.append(s.join([sensdate, senstime]))
        return sensdate_list

    def get_procdate(self):
        # returns the processing date from the manifest file
        root = et.fromstring(self.manifest)
        for elem in root.iter('metadataObject'):
            if elem.attrib['ID'] == 'processing':
                for subel in elem[0][0]:
                    dict_att = dict(subel.attrib)
        s = (' ')
        p_start = s.join([dict_att['start'][:10], dict_att['start'][11:]])
        return p_start

    def get_procvers(self):
         # returns the processor version from the manifest file
        root = et.fromstring(self.manifest)
        for elem in root.iter('metadataObject'):
            if elem.attrib['ID'] == 'processing':
                for subel in elem[0][0][0][0]:
                    dict_att = dict(subel.attrib)
        return dict_att['version']

    def get_linesPerBurst(self):
        # returns a list of lines per burst in each annotation file
        lpb_list = []
        for item in self.annotation:
            root = et.fromstring(item)
            lpb_list.append(np.int(root.find('swathTiming').find('linesPerBurst').text))
        return lpb_list    

    def get_pixelsPerBurst(self):
        # returns a list of pixels per burst in each annotation file
        ppb_list = []
        for item in self.annotation:
            root = et.fromstring(item)
            ppb_list.append(np.int(root.find('swathTiming').find('samplesPerBurst').text))
        return ppb_list

    def get_burstlist(self):
        # Returns the full xml element for the burstlist section for each of the annotation files
        bursts_list = []
        for item in self.annotation:
            root = et.fromstring(item)
            bursts_list.append(root.find('swathTiming').find('burstList'))
        return bursts_list

    def get_geolocGrid(self):
        # Returns the full xml element for the geolocation grid for each annotation file 
        glgrid_list = []
        for item in self.annotation:
            root = et.fromstring(item)
            glgrid_list.append(root.find('geolocationGrid')[0])
        return glgrid_list       

    def get_burstGeo_coords(self):
        # Returns two lists, one of the corners, one of the centres
        # Get the list of GeolocationGrids from each annotation file
        glggrid_list = self.get_geolocGrid()
        # Initialize list to store geocoordinates
        corners_list = []
        centres_list = []
        # Extract the edge coordinates of each burst
        for coordlist in glggrid_list:
            npts = coordlist.get('count')
                #print "Number of Geotagged points: " + npts
            pts = 0
            index = 0
            lo = [ ];  la = [ ]; ind = [ ]
            for pt in coordlist:
                pts = pts+1
                lon = pt.find('longitude').text
                lat = pt.find('latitude').text
                line = pt.find('line').text
                pixel = pt.find('pixel').text
                # The following assumes there are 21 points per line of geotagged points in xml files
                # If this breaks divide npts/(nbursts+1)
                if pts % 21 == 0 or pixel == "0":
                    index = index+1
                    lo = lo + [float(lon)]
                    la = la + [float(lat)]
                    ind = ind + [float(index)]

            #print lo, la, ind
            #print
            # Now have 3 lists which contain data for each relevant point in the coordinate list 
            # The data is ordered like [pixel0(corner1burst1), mod21(corner2burst1), pixel0(corner3burst1-corner1burst2), mod21(corner4burst1-corner2burst2), etc]
            # Now get the corners from these and make a list of corners for each one of the annotation files
            burstID = 0
            burstcorners_list = []
            burstcentres_list = []
            for i in range(3,len(lo)+1,2):
                corners = []
                burstID = burstID+1
                corners.append((la[i-3],lo[i-3])) # corner1, lat lon tuple
                corners.append((la[i-2],lo[i-2])) # corner2
                corners.append((la[i-1],lo[i-1])) # corner3
                corners.append((la[i],lo[i])) # corner4  
                burstcorners_list.append(corners) # Takes the list of corner coordinates and appends to a list of corners for all the bursts
                lat_list_deg = [la[i-3], la[i-2], la[i-1], la[i]]
                lon_list_deg = [lo[i-3], lo[i-2], lo[i-1], lo[i]]
                # Calculate central coordinates of burst
                # First convert the lat lons into radians
                lat_list_rad = [radians(x) for x in lat_list_deg]
                lon_list_rad = [radians(x) for x in lon_list_deg]
                # Convert the radians to cartesians (unit radius)
                x_list = []
                y_list = []
                z_list = []
                for j in range(len(lat_list_rad)):
                    x_list.append(cos(lat_list_rad[j]) * cos(lon_list_rad[j]))
                    y_list.append(cos(lat_list_rad[j]) * sin(lon_list_rad[j]))
                    z_list.append(sin(lat_list_rad[j]))
                # The average is not weighted and we assume we always have 4 vertices
                tot_weight = 4
                ave_x = fsum(x_list)/tot_weight
                ave_y = fsum(y_list)/tot_weight
                ave_z = fsum(z_list)/tot_weight
                # Convert the average cartesian coords back to radians
                ave_lon = atan2(ave_y, ave_x)
                hyp = sqrt(ave_x * ave_x + ave_y * ave_y)
                ave_lat = atan2(ave_z, hyp)
                # Convert back to degrees
                ave_lat_deg = degrees(ave_lat)
                ave_lon_deg = degrees(ave_lon)
                burstcentres_list.append((ave_lat_deg, ave_lon_deg)) # Takes a lat,lon tuple for the centre and appends to a list of centres for all the bursts
            #print burstcentres_list
            #print
            # Now add the list of burst corners and centres for the last annotation file onto the rest
            corners_list.append(burstcorners_list) 
            centres_list.append(burstcentres_list)
            #print centres_list
            #print
        return corners_list, centres_list

    def get_burst_ANXID(self):
        # returns a list of the ANX ID of each burst for each ann file
        burstlist_list = self.get_burstlist()
        annANX_list = []
        for burstlist in burstlist_list:
            burstANX_list = []
            for burstno, burstinfo in enumerate(burstlist):
                aziAnxTime = np.float32(burstinfo.find('azimuthAnxTime').text)
                burstid = np.int32(np.round(aziAnxTime*10))
                burstANX_list.append(burstid)
            annANX_list.append(burstANX_list)
        return annANX_list
                
    def get_mgrsID(self):
        # Returns a list of the MGRS ID of each burst for each annotation file
        corners_list, centres_list = self.get_burstGeo_coords()
        mgrs_ann_list = []
        for burstcentres_list in centres_list:
            #print burstcentres_list
            print()
            mgrs_burst_list = []
            for centre in burstcentres_list:
                m = mgrs.MGRS()
                mgrs_burst_list.append(m.toMGRS(centre[0], centre[1]))
            mgrs_ann_list.append(mgrs_burst_list)
        return mgrs_ann_list
    
