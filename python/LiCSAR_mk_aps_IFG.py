#!/usr/bin/env python

import sys
import numpy as np

# Load input variable date master and date slave
indir = sys.argv[1]
outdir = sys.argv[2]
m_date = sys.argv[3]
s_date = sys.argv[4]

#print '------------------------------------------------'
#print 'Load (read from file) delays maps from master and slave          '
#print '------------------------------------------------'
minfile = open(indir+'/'+m_date+'/'+m_date+'.aps','rb');
master  = np.fromfile( file=minfile , dtype=np.float32 , count=-1)
master.byteswap(True) # Because files are prepared for GAMMA, we have to swap_bytes
sinfile = open(indir+'/'+s_date+'/'+s_date+'.aps','rb');
slave   = np.fromfile( file=sinfile , dtype=np.float32 , count=-1)
slave.byteswap(True) # Because files are prepared for GAMMA, we have to swap_bytes
minfile.close()
sinfile.close()

#print '------------------------------------------------'
#print 'Compute delays_maps = slave - master            '
#print '------------------------------------------------'
diffdelay = slave - master
diffdelay = diffdelay - np.mean(diffdelay)

#print '------------------------------------------------'
#print 'Save to file: delays maps master_slave (aps.delay.map_4l.raw)       '
#print '------------------------------------------------'
outfile = open(outdir+'/'+m_date+'_'+s_date+'.aps','wb')
var2write = diffdelay.astype(np.float32)
var2write.byteswap(True) # Because files are prepared for GAMMA, we have to swap_bytes
var2write.tofile(outfile)
outfile.close()

#print 'APS maps for interferogram completed OK'
# END of the PROGRAM
#print '------------------------------------------------'



