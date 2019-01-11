#!/usr/bin/env python

import sys, zipfile, re
import xml.etree.ElementTree as et
import os.path

def usage():
  print """
### Read a zip file and return number of existing bursts per subswath

zipSAFE_2_NBursts.py
  Program reads a compressed SAFE/ file to estimate central coordinates of each burst
  
  Usage: zipSAFE_2_BurstList.py zipfile
    zipfile   (input) name of the compressed Sentinel-1 SAFE/ file

 Author: Pablo J. Gonzalez, [p.j.gonzalez@leeds.ac.uk]
Version: 1.0
   Date: 02-Nov-2015
"""

if len(sys.argv) < 2:
  usage()
  sys.exit(-1)
  
# Parse input 
zipfilepath=sys.argv[1]
if len(sys.argv) == 3:
  IW=sys.argv[2]
zfile=zipfile.ZipFile(zipfilepath, 'r') # Open the zip file

filebasename=os.path.splitext(os.path.basename(zipfilepath))[0]
datestr=filebasename[17:32]

# Identify the annotation files containing the coordinates of the bursts

if len(sys.argv) == 3:
  list_annfile=['annotation/s1.-iw'+IW+'-slc-vv']
else :
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
        
        # Initialize list to store geocoordinates
        lo = [ ];  la = [ ]; ind = [ ]
        # Extract the edge coordinates of each burst
        for swathTiming in root.findall('swathTiming'):
	  burst=swathTiming.find('burstList')
          print burst.get('count')

        # Catch exception if root not available
        if root.tag == "root":
          print "I'm the root"
        
      except et.ParseError:
        print(fnames)
        sys.exit(2)
     
zfile.close() #Close zipfile

  
      
