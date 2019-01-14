#!/usr/bin/env python

import sys

def usage():
  print("""
### Read input text file with series of lonmin lonmax latmin latmax coordinates and 
    output global lonmin lonmax latmin latmax coordinates v1.0 21-Oct-2015 PJG
  
SLCbb_minmax.py
  Program output coordinates of bounding box for long and latitude encompassing all bursts.
  
  Usage: SLCbb_minmax.py coordinates.txt
    coordinates.txt   (input) name of the text file with lonmin lonmax latmin latmax coordinates
        
""")

if len(sys.argv) < 2:
  usage()
  sys.exit(-1)

minlon = []
maxlon = []
minlat = []
maxlat = []
with open(sys.argv[1]) as f:
    for line in f:                   # loop over the rows
        fields = line.split()        # parse the columns
        rowdata = list(map(float, fields)) # convert text to numbers
         # This is hardcoded, but works for now.
        minlon = minlon + [rowdata[0]]
        maxlon = maxlon + [rowdata[1]]
        minlat = minlat + [rowdata[2]]
        maxlat = maxlat + [rowdata[3]]
        #data = data + [rowdata]         # accumulate the results
print(min(minlon), max(maxlon), min(minlat), max(maxlat))