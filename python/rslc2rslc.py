#!/usr/bin/env python

import sys
import LiCSAR_iofunc as LICSARio
#import subprocess
#import re 

def usage():
  print("""
### rslc2rslc.py
  Program reads a complex big endian file and an offset value
  to generate an output complex binary file in big endian of 
  the original file displaced offset number of pixels

  Usage: rslc2rslc.py file.rslc outfile.rslc offset [outtype]
    file.rslc   (input)  name of the big-endian complex float/short binary file
    outfile.rslc  (output) name of the big-endian output complex float/short binary file    
    offset        (input)  amount of pixels to displace the input file into the output file
    outtype    either FCOMPLEX or SCOMPLEX

rslc2rslc.py v1.0 14-Apr-2016 PJG
(2026: now used also to convert to another ?COMPLEX format
Part of LiCSAR software package
""")

if len(sys.argv) < 4:
  print(""" ERROR: Wrong number of input arguments """)
  usage()
  sys.exit(-1)

rslc=sys.argv[1]
outrslc=sys.argv[2]
azoffset=int(sys.argv[3])


width1  = int(LICSARio.width_slc_gamma_par(rslc))
length1 = int(LICSARio.length_slc_gamma_par(rslc))
dtype1  = LICSARio.dtype_slc_gamma_par(rslc)

try:
  dtype2=sys.argv[4]
except IndexError:
  dtype2 = dtype1
  print('output file will be in '+dtype2)

#print "Master image is: ",RSLCmasterfile, "(",width1,",",length1,")"
#print "Cropped slave image is: ",RSLCslavefile, "(",width2,",",length2,")"

LICSARio.shift_convert_fast_slc_gamma(
  rslc, outrslc, width1, length1, dtype1, dtype2, azoffset)
