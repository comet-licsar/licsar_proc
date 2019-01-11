#!/usr/bin/env python

import sys,os.path
import numpy as np
import LiCSAR_iofunc as LICSARio
#import subprocess
#import re 

def usage():
  print """
### rslc2rslc.py
  Program reads a complex big endian file and an offset value
  to generate an output complex binary file in big endian of 
  the original file displaced offset number of pixels

  Usage: slc2rslc.py rslc.file newrslc.file offset
    height.file   (input)  name of the big-endian complex float/short binary file
    geomask.file  (output) name of the big-endian output complex float/short binary file    
    offset        (input)  amount of pixels to displace the input file into the output file

rslc2rslc.py v1.0 14-Apr-2016 PJG
Part of LiCSAR software package
"""

if len(sys.argv) < 3:
  print """ ERROR: Wrong number of input arguments """
  usage()
  sys.exit(-1)

RSLCmasterfile=sys.argv[1]
RSLCslavefile=sys.argv[2]
azoffset=int(sys.argv[3])

width1  = int(LICSARio.width_slc_gamma_par(RSLCmasterfile))
length1 = int(LICSARio.length_slc_gamma_par(RSLCmasterfile))
dtype1  = LICSARio.dtype_slc_gamma_par(RSLCmasterfile)
dtype2  = LICSARio.dtype_slc_gamma_par(RSLCslavefile)

#print "Master image is: ",RSLCmasterfile, "(",width1,",",length1,")"
#print "Cropped slave image is: ",RSLCslavefile, "(",width2,",",length2,")"

# Read the slave image to be moved into the master geometry
if (dtype2 == 'SCOMPLEX'):
  data = LICSARio.read_fast_slc_gamma_scomplex(RSLCslavefile)
else:
  print " rslc2rslc.py encountered during reading an unsupported data type: ", dtype2
  #print "Size Slave of complex file after reading: ", data.shape
  #print "Data type/Numpy type: ", type(data), data.dtype

# Shift the slave image into the position in the size of the master image
outdata=np.zeros((length1,width1),np.complex64)
outdata[0+azoffset:data.shape[0]+azoffset,:]=data

# Write the result to disk
if (dtype2 == 'SCOMPLEX'):
  LICSARio.write_fast_slc_gamma_scomplex(outdata,RSLCslavefile[:-4]+'full.rslc')
else:
  print " rslc2rslc.py encountered during writing an unsupported data type: ", dtype2
