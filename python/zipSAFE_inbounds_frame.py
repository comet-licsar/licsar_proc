#!/usr/bin/env python

import sys, zipfile, re
import xml.etree.ElementTree as et

# Check if a point is inside a polygon
def point_in_poly(x,y,poly):
  n = len(poly)
  inside = False
  p1x,p1y = poly[0]
  for i in range(n+1):
    p2x,p2y = poly[i % n]
    if y > min(p1y,p2y):
      if y <= max(p1y,p2y):
        if x <= max(p1x,p2x):
          if p1y != p2y:
            xints = (y-p1y)*(p2x-p1x)/(p2y-p1y)+p1x
          if p1x == p2x or x <= xints:
            inside = not inside
    p1x,p1y = p2x,p2y
  return inside

def points2polygon(west,east,north,south):
  xy1= (west,north) # upper left
  xy2= (east,north) # upper right
  xy3= (east,south) # lower right
  xy4= (west,south) # lower left
  xy5= (west,north) # upper left
  polygon = [ xy1, xy2, xy3, xy4, xy5 ]
  return polygon

# Subroutine to output central coordinates from a SAFE zip_file
def pyC_zip2slc(zipfilepath):
  # Subroutine to output bounding box coordinates from each subswath of a SAFE zip_file
  # Output is a structure with three fields for subswaths and coordinates in it
  # slc.subswath[1,2,3].[ lon[1,...,Nburst], lat[1,...,Nburst] ]
  
  # Open the zip file
  zfile=zipfile.ZipFile(zipfilepath, 'r') 
  # Initialize slc structure
  lonmin=[]; lonmax=[]; latmin=[]; latmax=[];
  lo = [ ];  la = [ ];
  #slc = [[ ] for x in range(3)] # Allocate a list of 3 lists

  # Identify the annotation files containing the coordinates of the bursts
  list_annfile=['annotation/s1.-iw1-slc-vv', 'annotation/s1.-iw2-slc-vv', 'annotation/s1.-iw3-slc-vv']
  for fnames in zfile.namelist(): # Loop over all fileanames within the zip file
    # Reset the lo and la lists 
    lo=[];  la=[];
    for annfile in list_annfile: # Loop through list_annfile filenames
      if re.search(annfile,fnames): # If a zip filename matches one for the list_annfile, do...
        # Open the file and Read the xml structure
        try:
          doc = et.parse(zfile.open(fnames))
          root = doc.getroot()
          
          for geoGridinfo in root.findall('geolocationGrid'):
            for coordlist in geoGridinfo:
              npts = coordlist.get('count')
              pts = 0; #ind = 0;
              for pt in coordlist:
                pts = pts+1
                lon = pt.find('longitude').text
                lat = pt.find('latitude').text
                line = pt.find('line').text
                pixel = pt.find('pixel').text
                lo = lo + [float(lon)] # Add the new line lon value to the lo list
                la = la + [float(lat)] # Add the new line lat value to the la list
          
          # Extract the min,max long and lat for each subswath
          lonmin = lonmin + [min(lo)]
          lonmax = lonmax + [max(lo)]
          latmin = latmin + [min(la)]
          latmax = latmax + [max(la)]
          
        except et.ParseError:
          print(fnames)
          sys.exit(2)
  zfile.close() #Close zipfile
  return lonmin, lonmax, latmin, latmax

def usage():
  print """
### Read a zip file and extract lonmin/lonmax/latmin/latmax coordinates v1.0 21-Oct-2015 PJG'  
  
zipSAFE_2_SLCbb.py
  Program output coordinates of bounding box for long and latitude encompassing all bursts.
  
  Usage: zipSAFE_2_SLCbb.py SAFE_file_zip
    zipfile   (input) name of the compressed Sentinel-1 SAFE/ file
        
"""


if len(sys.argv) < 2:
  usage()
  sys.exit(-1)

# Parse input 
lonmin, lonmax, latmin, latmax = pyC_zip2slc(sys.argv[1])
west=min(lonmin)
east=max(lonmax) 
south=min(latmin) 
north=max(latmax)
zipedgepts=points2polygon(west,east,north,south)

frame_polygon = []
with open(sys.argv[2]) as f:
  for line in f:
    frame_polygon.append([float(n) for n in line.strip().split(' ')])

numpts_within_polygon=0
# Loop over bursts in master
for index in range(0,len(zipedgepts),1):
    if ( point_in_poly(zipedgepts[index][0],zipedgepts[index][1],frame_polygon) ):
        numpts_within_polygon=numpts_within_polygon+1

if numpts_within_polygon > 0:
    print "1"
else:
    print "0"


