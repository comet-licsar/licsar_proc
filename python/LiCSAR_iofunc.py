#!/usr/bin/env python

import numpy as np
import subprocess

########################################################################
########################################################################
# Read metadata from gamma parameter files
######################################################################## 
def width_slc_gamma_par(file):
  calling=['grep range_samples: '+file+'.par | awk \'{print $2}\' ']
  child=subprocess.Popen(calling,stdout=subprocess.PIPE,shell=True)
  width=child.communicate()[0]
  return width[:-1]
def length_slc_gamma_par(file):
  calling=['grep azimuth_lines: '+file+'.par | awk \'{print $2}\' ']
  child=subprocess.Popen(calling,stdout=subprocess.PIPE,shell=True)
  length=child.communicate()[0]
  return length[:-1]
def dtype_slc_gamma_par(file):
  calling=['grep image_format: '+file+'.par | awk \'{print $2}\' ']
  child=subprocess.Popen(calling,stdout=subprocess.PIPE,shell=True)
  dtype=child.communicate()[0]
  return dtype[:-1]  
########################################################################
########################################################################

########################################################################
########################################################################
# Read/Write Float type
########################################################################
def read_slc_gamma_float(file):
  width  = int(width_slc_gamma_par(file))
  length = int(length_slc_gamma_par(file))
  data = np.fromfile(file,np.float32,length*width).reshape(length,width)
  return data
########################################################################
def write_slc_gamma_float(data,outname):
   data=np.array(data,dtype=np.float32)
   data.tofile(outname)  
########################################################################
########################################################################

########################################################################
########################################################################
# Read/Write Complex int16 type
########################################################################
def read_slc_gamma_scomplex(file):
  # read gamma scomplex data, i.e. .slc file.  
  wid1 = int(width_slc_gamma_par(file))
  len2 = int(length_slc_gamma_par(file))
  data = np.fromfile(file,np.int16,len2*2*wid1)
  id1 = range(0,2*len2*wid1,2)
  id2 = range(1,2*len2*wid1,2)
  real = data[id1].reshape(len2,wid1)
  imag = data[id2].reshape(len2,wid1)
  data_cpx = real + imag*1j
  return data_cpx
def read_fast_slc_gamma_scomplex(file):
  # read faster a gamma scomplex data, i.e. .slc file.    
  wid1 = int(width_slc_gamma_par(file))
  len2 = int(length_slc_gamma_par(file))
  data = np.fromfile(file,np.int16,len2*2*wid1)
  real   = data[::2].reshape(len2,wid1)
  imag   = data[1::2].reshape(len2,wid1)
  data_cpx = np.zeros(real.shape,dtype=np.complex64)
  data_cpx.real = real
  data_cpx.imag = imag
  return data_cpx
########################################################################
def write_slc_gamma_scomplex(data,outname):
   # write gamma scomplex data, i.e. .slc file.
   nlines = data.shape[0]
   WIDTH  = data.shape[1]
   id1 = range(0,2*nlines*WIDTH,2)
   id2 = range(1,2*nlines*WIDTH,2)
   F=np.zeros([2*nlines*WIDTH,1],np.int16)
   F[id1]=np.reshape(data.real,(nlines*WIDTH,1))
   F[id2]=np.reshape(data.imag,(nlines*WIDTH,1))
   F.tofile(outname)
def write_fast_slc_gamma_scomplex(data,outname):
   # write gamma scomplex data, i.e. .slc file.
   nlines = data.shape[0]
   WIDTH  = data.shape[1]
   F=np.zeros([2*nlines*WIDTH,1],np.int16)
   F[::2]=np.reshape(data.real,(nlines*WIDTH,1))
   F[1::2]=np.reshape(data.imag,(nlines*WIDTH,1))
   F.tofile(outname)  
########################################################################
########################################################################

########################################################################
########################################################################
# Read/Write Complex float32 type
########################################################################
def read_fast_slc_gamma_fcomplex(file):
  # read faster a gamma scomplex data, i.e. .slc file.    
  wid1 = int(width_slc_gamma_par(file))
  len2 = int(length_slc_gamma_par(file))
  data = np.fromfile(file,np.float32,len2*2*wid1)
  real   = data[::2].reshape(len2,wid1)
  imag   = data[1::2].reshape(len2,wid1)
  data_cpx = np.zeros(real.shape,dtype=np.complex64)
  data_cpx.real = real
  data_cpx.imag = imag
  return data_cpx  
def read_fcomplex(file,wid1,len2):
  # read fcomplex data, i.e. .diff file.    
  data_cpx = np.fromfile(file,np.complex64,len2*2*wid1).reshape(len2,wid1)
  return data_cpx   
########################################################################
def write_fast_slc_gamma_fcomplex(data,outname):
   # write gamma scomplex data, i.e. .slc file.
   nlines = data.shape[0]
   WIDTH  = data.shape[1]
   F=np.zeros([2*nlines*WIDTH,1],np.float32)
   F[::2]=np.reshape(data.real,(nlines*WIDTH,1))
   F[1::2]=np.reshape(data.imag,(nlines*WIDTH,1))
   F.tofile(outname)
def write_fcomplex(data,outname):
   # write fcomplex data, i.e. .diff file.
   nlines = data.shape[0]
   WIDTH  = data.shape[1]
   F=np.zeros([2*nlines*WIDTH,1],np.float32)
   F[::2]=np.reshape(data.real,(nlines*WIDTH,1))
   F[1::2]=np.reshape(data.imag,(nlines*WIDTH,1))
   F.tofile(outname)     
########################################################################
########################################################################



def realimag2cpxint16(real,imag):
  ci2type = np.dtype([('re', np.int16), ('im', np.int16)])
  data_cpx = np.zeros(real.shape,ci2type)
  data_cpx['re'] = real 
  data_cpx['im'] = imag #imag*1j
  return data_cpx


 