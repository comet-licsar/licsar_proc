#!/usr/bin/env python

import sys,os.path
import numpy as np
import LiCSAR_iofunc as LICSARio

def usage():
  print """
### IFGphaseoffset.py
  Program reads a complex big endian file and add a phase offset value, X
  
  Optionally, if input a mosaicked SLC, one can apply the phase offset (X)
  to only the left portion of the image until the column "colIW1"
                         [XXXXXXXX|---------------]
                         [XXXXXXXX|---------------]
                         [XXXXXXXX|---------------]
                         [---colIW1----------width]    

  Usage: IFGphaseoffset.py input.ifg width length output.ifg phshift [colIW1]"
    input.slc   (input)  name of the big-endian complex float binary file
    width       (input)  width (num of columns) in complex float binary file
    length      (input)  length (num of lines) in complex float binary file
    output.slc  (output) name of the big-endian output complex float/short binary file    
    offset      (input)  phase offset (in radians)
    colIW1      (input)  [optional] max. column to apply phase offset [meant for mosaics]

IFGphaseoffset.py v1.0 28-Apr-2016 PJG
Part of LiCSAR software package
"""

if (len(sys.argv) < 6) or (len(sys.argv) > 7) :
  print """ ERROR: Wrong number of input arguments """
  usage()
  sys.exit(-1)

IFGinfile=sys.argv[1]
width1=int(sys.argv[2])
length1=int(sys.argv[3])
IFGoutfile=sys.argv[4]
phoffset=float(sys.argv[5])
if len(sys.argv) == 7:
  colIW1=int(sys.argv[6])
else:
  colIW1=width1

cpxIFG = LICSARio.read_fcomplex(IFGinfile,width1,length1)
cpxIFG.byteswap(True) # Swap bytes - from big endian to little endian

# Shift phase of the SLC 
cph = np.cos(phoffset)
sph = np.sin(phoffset)
real1 = cpxIFG[:,0:colIW1].real * cph - cpxIFG[:,0:colIW1].imag * sph
imag1 = cpxIFG[:,0:colIW1].real * sph + cpxIFG[:,0:colIW1].imag * cph
cpxIFG[:,0:colIW1] = real1 + imag1*1j

cpxIFG.byteswap(True) # Swap bytes - from big endian to little endian
LICSARio.write_fcomplex(cpxIFG,IFGoutfile)






