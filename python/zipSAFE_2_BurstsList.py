#!/usr/bin/env python

import sys, zipfile, re
import xml.etree.ElementTree as et
import os.path

def usage():
  print("""
### Read a zip file and extract info for existing bursts

zipSAFE_2_BurstList.py
  Program reads a compressed SAFE/ file to estimate central coordinates of each burst
  
  Usage: zipSAFE_2_BurstList.py zipfile
    zipfile   (input) name of the compressed Sentinel-1 SAFE/ file

 Author: Pablo J. Gonzalez, [p.j.gonzalez@leeds.ac.uk]
Version: 1.0
   Date: 02-Nov-2015
""")

if len(sys.argv) < 2:
  usage()
  sys.exit(-1)
  
# Parse input 
zipfilepath=sys.argv[1]
zfile=zipfile.ZipFile(zipfilepath, 'r') # Open the zip file

filebasename=os.path.splitext(os.path.basename(zipfilepath))[0]
datestr=filebasename[17:32]

# Identify the annotation files containing the coordinates of the bursts
list_annfile=['annotation/s1.-iw1-slc-vv', 'annotation/s1.-iw2-slc-vv', 'annotation/s1.-iw3-slc-vv']
for fnames in zfile.namelist():
  for annfile in list_annfile:
    if re.search(annfile,fnames):
      # Open the file and Read the xml structure
      try:
        doc = et.parse(zfile.open(fnames))
        root = doc.getroot()
        # Extract the subswath info
        for headerinfo in root.findall('adsHeader'):
          swath=headerinfo.find('swath').text
          #relativeOrbit=headerinfo.find('absoluteOrbitNumber').text
        #file0 = open("T"+relativeOrbit+"_IW"+swath+".burstlist", "w")
        file0 = open(datestr+"_"+swath+".burstlist", "w")
        
        # Initialize list to store geocoordinates
        lo = [ ];  la = [ ]; ind = [ ]
        # Extract the edge coordinates of each burst
        for geoGridinfo in root.findall('geolocationGrid'):
          for coordlist in geoGridinfo:
            npts = coordlist.get('count')
            #print "Number of Geotagged points: " + npts
            pts = 0
            index = 0
            for pt in coordlist:
              pts = pts+1
              lon = pt.find('longitude').text
              lat = pt.find('latitude').text
              line = pt.find('line').text
              pixel = pt.find('pixel').text
              # The following assumes there are 21 points per line of geotagged points in xml files
              if pts % 21 == 0 or pixel == "0":
                index = index+1
                lo = lo + [float(lon)]
                la = la + [float(lat)]
                ind = ind + [float(index)]
                
        burstID = 0
        for i in range(3,len(lo)+1,2):
          burstID = burstID+1
          # Calculate central coordinates of burst
          lon0 = (lo[i-3]+lo[i-2]+lo[i-1]+lo[i])/4
          lat0 = (la[i-3]+la[i-2]+la[i-1]+la[i])/4
          # Write central coordinates
          file0.write("%s %s %s\n" % (lon0, lat0, burstID) )

        # Close files
        file0.close()
        
        # Catch exception if root not available
        if root.tag == "root":
          print("I'm the root")
        # Close file where we put the coordinate list
        #text_file.close()
        
      except et.ParseError:
        print(fnames)
        sys.exit(2)
zfile.close() #Close zipfile

  
      
