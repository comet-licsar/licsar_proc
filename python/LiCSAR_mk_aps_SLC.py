#!/usr/bin/env python

import sys
import numpy as np
import pyaps as pa

#LiCSAR_mk_aps_SLC.py $slcdir $outdir $slcdate geo/EQA.dem.littleendian $hour $lambda geo/${masterdate}.deg.inc

# Load input variable hour - closest day hour time to SAR acquisition time
if len(sys.argv) == 6: 
  outdir  = sys.argv[1]
  slcdate = sys.argv[2]
  demfilename = sys.argv[3]
  hour    = sys.argv[4]
  radarwv = float(sys.argv[5])
  incfile = 23.0
else:  
  outdir  = sys.argv[1]
  slcdate = sys.argv[2]
  demfilename = sys.argv[3]
  hour    = sys.argv[4]
  radarwv = float(sys.argv[5])
  incfile = sys.argv[6]
  
pa.ECMWFdload([str(slcdate)],str(hour),str(outdir))
aps1 = pa.PyAPS_geo(str(outdir)+'/'+'ERA-Int_'+str(slcdate)+'_'+str(hour)+'.grb',demfilename,grib='ECMWF',verb=True,demfmt='HGT',demtype=np.float32)
aps1.getdelay(outdir+'/'+str(slcdate)+'.geo.delay',inc=incfile,wvl=radarwv)




