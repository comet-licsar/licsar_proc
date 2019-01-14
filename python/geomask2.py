#!/usr/bin/env python

import sys,os.path
import numpy as np

def usage():
  print("""
### geomask2.py
  Program reads two big endian files with elevation and overlay/shadow pixels 
  to generate mask output a binary file in big endian

  Usage: geomask2.py height.file ls_map.file width lenght 
    height.file   (input)  name of the big-endian elevation height float binary file
    ls_map.file   (input)  name of the big-endian layover and shadow map float binary file    
    geomask.file  (output) name of the big-endian output mask float binary file    
    width         (input)  width of the binary matrix
    length        (input)  length of the binary matrix

geomask2.py v1.0 22-Mar-2016 PJG
Part of LiCSAR software package
""")

if len(sys.argv) < 5:
  print(""" ERROR: Wrong number of input arguments """)
  usage()
  sys.exit(-1)

heightfile=sys.argv[1]
shadowfile=sys.argv[2]
outfile=sys.argv[3]
nx=int(sys.argv[4]) # Width of the files [number of columns]
ny=int(sys.argv[5]) # length of the files [number of rows]

# Define input data depth
filetype=np.float32
filetype2=np.int8

assert os.path.isfile(heightfile), 'LiCSAR - geomask.py: height (elevations) angle file does not exist'
assert os.path.isfile(shadowfile), 'LiCSAR - geomask.py: layover/shadow pixels file does not exist'
  
fin = open(heightfile,'rb'); height = np.fromfile(file=fin,dtype=filetype,count=nx*ny).reshape(ny,nx); fin.close()
height.byteswap(True) # Swap bytes - from big endian to little endian

fin = open(shadowfile,'rb'); shadow = np.fromfile(file=fin,dtype=filetype2,count=nx*ny).reshape(ny,nx); fin.close()
shadow.astype(np.float32)
shadow_masked=np.where(shadow>1., 0., 1.)  # Mask any value larger than 1, if condition is true 0 else 1

out = np.multiply(shadow_masked, height); 
out = out.astype(np.float32)
out1= out.byteswap(True); 

outFile = isinstance(outfile,str)
if outFile:
  fout = open(outfile,'wb')
out1.tofile(fout)
fout.close()
