#!/usr/bin/env python

import sys, math

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

def commonbursts(master,polygon):
  # Compute distance between central coordinates for each burst to figure out if the burst are the same
  # If distance is larger than (~5 km), I assume they are different bursts. 
  # Interburst distance is crudely approximate for latitude and more precisely computed for longitudes
  
  # Initialize list of lists to store indexes of common bursts
  bn = [ ]; # Burst number
  # Loop over bursts in master
  for m_burstID in range(0,len(master),1):
    if ( point_in_poly(master[m_burstID][0],master[m_burstID][1],polygon) ): # if ( point_in_poly(master[m_burstID][0],master[m_burstID][1],polygon) ):
      print(master[m_burstID][0] , master[m_burstID][1])
  
def usage():
  print("""
### Read one burst coordinate lists and filter by burst centre coordinates within a given polygon, 
    reporting bursts list coordinates. 
  
    BurstListCoordInPolyFile.py
    Program compares coordinates of burst centers filtered by the ones falling within a given polygon file.
    
    Usage: BurstListCoordInPolyFile.py master_burstlist polyfile
      master_burstlist  (input) filename master burst list with format [long, lat, burstID]
              polyfile  (input) polyfile consisting of two columns of [long, lat], defining region. 

 Author: Pablo J. Gonzalez, [p.j.gonzalez@leeds.ac.uk]
Version: 1.0
   Date: 3-Dec-2015
""")


if len(sys.argv) < 3:
  usage()
  sys.exit(-1)

master = []
with open(sys.argv[1]) as f:
  for line in f:
    master.append([float(n) for n in line.strip().split(' ')])

polygon = []
with open(sys.argv[2]) as f:
  for line in f:
    polygon.append([float(n) for n in line.strip().split(' ')])

commonbursts(master,polygon)


