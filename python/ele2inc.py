#!/usr/bin/env python

import sys,os.path
import numpy as np

def usage():
  print("""
### Convert from elevation angles in radians to incidence angles in degrees v1.0 28-Jan-2016 PJG'  
  
ele2inc.py
  Program reads a big endian file with output a binary file in big endian
  
  Usage: ele2inc.py elev.file inc.file
    elev.file   (input)  name of the big-endian elevation angles float binary file
    inc.file    (output) name of the big-endian incidence angles float binary file
    width       (input)  width of the binary matrix
    length      (input)  length of the binary matrix
    
""")

if len(sys.argv) < 4:
  usage()
  sys.exit(-1)

elevfile=sys.argv[1]
incfile=sys.argv[2]
nx=int(sys.argv[3]) # Width of the files [number of columns]
ny=int(sys.argv[4]) # length of the files [number of rows]

# Define input data depth
filetype=np.float32

assert os.path.isfile(elevfile), 'LiCSAR - ele2inc.py: elevation angle file does not exist'
hfile = elevfile
  
fin = open(hfile,'rb');
elevrad = np.fromfile(file=fin,dtype=filetype,count=nx*ny).reshape(ny,nx)
elevrad.byteswap(True) # Swap bytes - from big endian to little endian
fin.close()

incdeg  = 90.0 - (elevrad*(180./np.pi)) # 90 - (rad2deg)
incdeg.byteswap(True)  # Swap bytes - from little endian to big endian

outFile = isinstance(incfile,str)
if outFile:
  fout = open(incfile,'wb')

incdeg.tofile(fout)
fout.close()

#dout = res.astype(np.float32)
#dout.tofile(fout)
#fout.close()

