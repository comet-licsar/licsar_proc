#!/usr/bin/env python

import numpy as np
import subprocess
from LiCSAR_lib.LiCSAR_misc import *

########################################################################
########################################################################
# Read metadata from gamma parameter files
######################################################################## 
def width_slc_gamma_par(file):
  res = grep1('range_samples:',file+'.par')
  res = res.split(':')[1].strip()
  return res
def length_slc_gamma_par(file):
  res = grep1('azimuth_lines: ',file+'.par')
  res = res.split(':')[1].strip()
  return res
def dtype_slc_gamma_par(file):
  res = grep1('image_format: ',file+'.par')
  res = res.split(':')[1].strip()
  return res
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
  id1 = list(range(0,2*len2*wid1,2))
  id2 = list(range(1,2*len2*wid1,2))
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
   id1 = list(range(0,2*nlines*WIDTH,2))
   id2 = list(range(1,2*nlines*WIDTH,2))
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

def shift_convert_fast_slc_gamma(infile, outname, width, length,
                                 input_type, output_type, azimuth_offset=0,
                                 chunk_lines=256):
  """Shift and convert a GAMMA complex raster using bounded row blocks."""
  import os

  types = {
      'SCOMPLEX': np.dtype(np.int16),
      'FCOMPLEX': np.dtype(np.float32),
  }
  if input_type not in types:
    raise ValueError('Unsupported input complex type: {}'.format(input_type))
  if output_type not in types:
    raise ValueError('Unsupported output complex type: {}'.format(output_type))
  if width <= 0 or length <= 0:
    raise ValueError('Raster width and length must be positive')
  if chunk_lines <= 0:
    raise ValueError('chunk_lines must be positive')
  if os.path.abspath(infile) == os.path.abspath(outname):
    raise ValueError('Input and output files must be different')

  input_dtype = types[input_type]
  output_dtype = types[output_type]
  row_values = 2 * width

  with open(infile, 'rb') as source, open(outname, 'wb') as destination:
    for dst_start in range(0, length, chunk_lines):
      dst_end = min(dst_start + chunk_lines, length)
      output = np.zeros((dst_end - dst_start, row_values), dtype=output_dtype)

      src_start = max(0, dst_start - azimuth_offset)
      src_end = min(length, dst_end - azimuth_offset)
      if src_start < src_end:
        source.seek(src_start * row_values * input_dtype.itemsize)
        values = np.fromfile(source, dtype=input_dtype,
                             count=(src_end - src_start) * row_values)
        expected = (src_end - src_start) * row_values
        if values.size != expected:
          raise IOError('Input complex raster is shorter than its parameter dimensions')
        dst_offset = src_start + azimuth_offset - dst_start
        output[dst_offset:dst_offset + (src_end - src_start)] = values.reshape(
            src_end - src_start, row_values).astype(output_dtype, copy=False)

      output.tofile(destination)
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


 
