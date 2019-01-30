#!/usr/bin/env python

import sys,os.path
import numpy as np
import LiCSAR_iofunc as LICSARio

def usage():
  print("""
### SLCphaseoffset.py
  Program reads a complex big endian file and add a phase offset value, X
  
  Optionally, if input a mosaicked SLC, one can apply the phase offset (X)
  to only the left portion of the image until the column "colIW1"
                         [XXXXXXXX|---------------]
                         [XXXXXXXX|---------------]
                         [XXXXXXXX|---------------]
                         [---colIW1----------width]    

  Usage: SLCphaseoffset.py input.slc output.slc phshift [colIW1]"
    input.slc   (input)  name of the big-endian complex float/short binary file
    output.slc  (output) name of the big-endian output complex float/short binary file    
    offset      (input)  phase offset (in radians)
    colIW1      (input)  [optional] max. column to apply phase offset [meant for mosaics]

SLCshiftphase.py v1.0 27-Apr-2016 PJG
Part of LiCSAR software package
""")

if len(sys.argv) < 3:
  print(""" ERROR: Wrong number of input arguments """)
  usage()
  sys.exit(-1)

SLCinfile=sys.argv[1]
SLCoutfile=sys.argv[2]
phoffset=float(sys.argv[3])
colIW1=int(sys.argv[4])

print("Input: ",SLCinfile,"  Output: ",SLCoutfile,"  Phase Offset: ", phoffset, " ColumnIW1: ", colIW1)
#if len(sys.argv) == 4:
#  colIW1=int(sys.argv[4])
#else:
#  colIW1=int(LICSARio.width_gamma_par(SLCinfile))

width1  = int(LICSARio.width_slc_gamma_par(SLCinfile))
length1 = int(LICSARio.length_slc_gamma_par(SLCinfile))
dtype1  = LICSARio.dtype_slc_gamma_par(SLCinfile)

#print "Master image is: ",RSLCmasterfile, "(",width1,",",length1,")"
#print "Cropped slave image is: ",RSLCslavefile, "(",width2,",",length2,")"

# Read the slave image to be moved into the master geometry
if (dtype1 == 'SCOMPLEX'):
  cpxSLC = LICSARio.read_fast_slc_gamma_scomplex(SLCinfile)
  cpxSLC.byteswap(True) # To operate with those matrices we have to swapbytes!!
else:
  print(" SLCshiftphase.py encountered during reading an unsupported data type: ", dtype1)

# Shift phase of the SLC 
cph = np.cos(phoffset)
sph = np.sin(phoffset)
real1 = cpxSLC[:,0:colIW1].real * cph - cpxSLC[:,0:colIW1].imag * sph
imag1 = cpxSLC[:,0:colIW1].real * sph + cpxSLC[:,0:colIW1].imag * cph
cpxSLC[:,0:colIW1].real = real1
cpxSLC[:,0:colIW1].imag = imag1
#cpxSLC[:,0:colIW1] = real1 + imag1*1j  

#cpxSLC[:,0:colIW1] = np.abs(cpxSLC[:,0:colIW1]) * np.exp(1j*(np.angle(cpxSLC[:,0:colIW1])+phoffset))

# Write the result to disk
if (dtype1 == 'SCOMPLEX'):
  cpxSLC.byteswap(True) 
  LICSARio.write_fast_slc_gamma_scomplex(cpxSLC,SLCoutfile)
else:
  print(" SLCphaseoffset.py encountered during writing an unsupported data type: ", dtype1)


#cpxSLC[:,0:colIW1] = np.abs(cpxSLC[:,0:colIW1])*np.exp(1j*(np.angle(cpxSLC[:,0:colIW1])*phaseoffset))
#cpxSLC[:,0:colIW1] = np.multiply(np.abs(cpxSLC[:,0:colIW1]),np.exp(1j*(np.angle(cpxSLC[:,0:colIW1])+phoffset)))

#c1 = cos(ph);
##s1 = -sin(ph);    /* conjugate phase to subtract from int-1 */
#s1 = sin(ph);
#diff[2*i] =   int1[2*i] * c1 - int1[2*i+1] * s1;
#diff[2*i+1] = int1[2*i] * s1 + int1[2*i+1] * c1;


