#!/usr/bin/env python

import sys, math

def rad2deg(radians):
  pi = math.pi
  degrees = 180 * radians / pi
  return degrees
  
def deg2rad(degrees):
  pi = math.pi
  radians = pi * degrees / 180
  return radians

#def points2polygon(west,east,north,south):
  #xy1= (west,north) # upper left
  #xy2= (east,north) # upper right
  #xy3= (east,south) # lower right
  #xy4= (west,south) # lower left
  #xy5= (west,north) # upper left
  #polygon = [ xy1, xy2, xy3, xy4, xy5 ]
  #return polygon

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

def commonbursts(master,slave,polygon):
  # Compute distance between central coordinates for each burst to figure out if the burst are the same
  # If distance is larger than (~5 km), I assume they are different bursts. 
  # Interburst distance is crudely approximate for latitude and more precisely computed for longitudes
  
  # Some constants to compute distances in km from comparing coordinates in decimal degrees
  a    = 6378.137 # equatorial radius of Earth [km]
  b    = 6356.752 # polar radius of Earth [km]
  dLat =  111.000; # Average value to map distance in latitude to kms [varies from 110.57km@0deg to 111.694km@90deg]
  ecc2 = ( pow(a,2) - pow(b,2) ) / pow(a,2) # Squared eccentricity
  # Initialize list of lists to store indexes of common bursts
  bn = [ ]; # Burst number
  # Loop over bursts in master
  for m_burstID in range(0,len(master),1):
    if ( point_in_poly(master[m_burstID][0],master[m_burstID][1],polygon) ): # if ( point_in_poly(master[m_burstID][0],master[m_burstID][1],polygon) ):
      # Loop over bursts in slave
      for s_burstID in range(0,len(slave),1): 
        if ( point_in_poly(slave[s_burstID][0],slave[s_burstID][1],polygon) ):
          phi1_rad = deg2rad(master[m_burstID][0])
          phi2_rad = deg2rad(slave[s_burstID][0])
          dx1 = (math.pi*a*math.cos(phi1_rad)) / (180*math.sqrt(1 - ecc2*pow(math.sin(phi1_rad),2)))
          dx2 = (math.pi*a*math.cos(phi2_rad)) / (180*math.sqrt(1 - ecc2*pow(math.sin(phi2_rad),2)))
          dx = pow((dx1-dx2), 2)
          dy = pow((master[m_burstID][1] - slave[s_burstID][1])*dLat, 2)
          d  = math.sqrt(dx+dy) # Distance in km
          if d<5: # If distance is smaller than 5 km I assume they belong to the same burst
            bn.append(m_burstID+1)
  if bn:
    bn1 = bn[0]  # 1st Burst ID for the master that matches the slave ones
    bn2 = bn[-1] # Last Burst ID for the master that matches the slave ones
    return bn1, bn2
  else:
    bn1 = 0  # 1st Burst ID for the master that matches the slave ones
    bn2 = 0  # Last Burst ID for the master that matches the slave ones
    return bn1, bn2

def usage():
  print """
### Read two burst coordinate lists and extract common bursts (filtered by burst 
centre coordinates within a given polygon), reporting first and last common bursts. 
v0.1 3-Dec-2015 PJG'  
  
BurstList2CommonGeoBurstsPolyFile.py
  Program compares coordinates of burst centers between a pair of S1 images 
  and output a burst list (numbers) with respect to burst order in the master.
  filtered by the ones falling within a given polygon file.
    
  Usage: BurstList2CommonGeoBurstsPolyFile.py master_burstlist slave_burstlist polyfile
    master_burstlist  (input) filename master burst list with format [long, lat, burstID]
     slave_burstlist  (input) filename slave burst list with format [long, lat, burstID]
            polyfile  (input) polyfile consisting of two columns of [long, lat], defining
                              region. 

 Author: Pablo J. Gonzalez, [p.j.gonzalez@leeds.ac.uk]
Version: 1.0
   Date: 3-Dec-2015
"""


if len(sys.argv) < 4:
  usage()
  sys.exit(-1)

master = []
with open(sys.argv[1]) as f:
  for line in f:
    master.append([float(n) for n in line.strip().split(' ')])
slave = []
with open(sys.argv[2]) as f:
  for line in f:
    slave.append([float(n) for n in line.strip().split(' ')])

polygon = []
with open(sys.argv[3]) as f:
  for line in f:
    polygon.append([float(n) for n in line.strip().split(' ')])

bn1,bn2 = commonbursts(master,slave,polygon)
print bn1, bn2 

