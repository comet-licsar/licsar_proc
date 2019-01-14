#!/usr/bin/env python

import sys,os.path
import numpy as np

def usage():
  print("""
### orient2heading.py
  Convert from orientation angles in radians to heading angles in degrees
  Program reads a big endian file with output a binary file in big endian
  
  Usage: orient2heading.py elev.file inc.file width length 
    orient.file   (input)  name of the big-endian orientation angles float binary file
    heading.file  (output) name of the big-endian heading angles float binary file
    width         (input)  width of the binary matrix
    length        (input)  length of the binary matrix

orient2heading.py v1.0 28-Jan-2016 PJG
Part of LiCSAR software package
""")

if len(sys.argv) < 4:
  usage()
  sys.exit(-1)

orientfile=sys.argv[1]
headingfile=sys.argv[2]
nx=int(sys.argv[3]) # Width of the files [number of columns]
ny=int(sys.argv[4]) # length of the files [number of rows]

# Define input data depth
filetype=np.float32

assert os.path.isfile(orientfile), 'LiCSAR - orient2heading.py: orientation angle file does not exist'
  
fin = open(orientfile,'rb');
orientrad = np.fromfile(file=fin,dtype=filetype,count=nx*ny).reshape(ny,nx)
orientrad.byteswap(True) # Swap bytes - from big endian to little endian
fin.close()

headingdeg  = -180.0 - (orientrad*(180./np.pi)) # 180 - (rad2deg)
headingdeg.byteswap(True)  # Swap bytes - from little endian to big endian

outFile = isinstance(headingfile,str)
if outFile:
  fout = open(headingfile,'wb')

headingdeg.tofile(fout)
fout.close()


