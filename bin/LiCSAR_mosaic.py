#!/usr/bin/env python
"""
    LiCSAR_mosaic - simple tool to remosaic a rslc.

    usage:
    LiCSAR_mosiac.py <procDir> <slaveDate> <masterDate> <rglks> <azlks>
"""
################################################################################
#imports
################################################################################
import sys
import datetime as dt
from LiCSAR_lib.coreg_lib import rebuild_rslc

################################################################################
#input args
################################################################################
procDir = sys.argv[1]
slaveDate = dt.datetime.strptime(sys.argv[2],'%Y%m%d')
masterDate = dt.datetime.strptime(sys.argv[3],'%Y%m%d')
rglks = int(sys.argv[4])
azlks = int(sys.argv[5])

################################################################################
#Call rebuild
################################################################################
rebuild_rslc(procDir,slaveDate,masterDate,rglks,azlks)
