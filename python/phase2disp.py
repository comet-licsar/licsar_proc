#!/usr/bin/env python

import sys,os.path
import numpy as np

def usage():
  print("""
### phase2disp.py
  Program reads a big endian file with unwrapped phases in radians 
  to generate a displacement map output to a binary file in big endian

  Usage: phase2disp.py unw.file lambda disp.file width lenght 
    unw.file   (input)  name of the big-endian unwrap phase float binary file [radians]
    lambda     (input)  radar wavelength [cm]
    disp.file  (output) name of the big-endian displacement map float binary file [cm]
    width      (input)  width of the binary matrix
    length     (input)  length of the binary matrix

phase2disp.py v1.0 10-Feb-2016 PJG
Part of LiCSAR software package
""")

if len(sys.argv) < 5:
  print(""" ERROR: Wrong number of input arguments """)
  usage()
  sys.exit(-1)

unwfile=sys.argv[1]
lam=float(sys.argv[2])
outfile=sys.argv[3]
nx=int(sys.argv[4]) # Width of the files [number of columns]
ny=int(sys.argv[5]) # length of the files [number of rows]

# Define input data depth
filetype=np.float32

assert os.path.isfile(unwfile), 'LiCSAR - geomask.py: height (elevations) angle file does not exist'
  
fin = open(unwfile,'rb');
unw = np.fromfile(file=fin,dtype=filetype,count=nx*ny).reshape(ny,nx)
unw.byteswap(True) # Swap bytes - from big endian to little endian
fin.close()

out = (lam/(4*np.pi))*unw       # ( (4 PI) / $lambda ) * variable
out.byteswap(True)  # Swap bytes - from little endian to big endian

outFile = isinstance(outfile,str)
if outFile:
  fout = open(outfile,'wb')

out.tofile(fout)
fout.close()

