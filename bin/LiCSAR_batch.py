#!/usr/bin/env python

import os
import subprocess as subp
import sys
import datetime as dt
try:
    import global_config as gc
except:
    print "Can't find local copy of global_config.py. You should copy it to your working directory."
    raise 
sys.path.insert(0,gc.LiCSAR_PATH+'/lib')
from LiCSAR_01_mk_images import cd
import LiCSAR_01_mk_images as L01mk
import LiCSAR_02_coreg as L02co
from LiCSAR_03_mk_ifgs import make_interferogram, make_baselines, make_interferograms_list
from LiCSAR_04_unwrap import main as unwrap_all

from LiCSAR_00_check_images import main as data_check

def main(argv=None):
    data_check()

    mdate = dt.date(int(gc.master[:4]),int(gc.master[4:6]),int(gc.master[6:8]))
    bperpfile = gc.procdir+'bperp_'+gc.master+'_file'
    
    mkargs = ['null', '-d', os.path.abspath(gc.procdir), '-y', gc.batchflag, '-p', gc.polygonfile, '-z', gc.ziplistfile, '-m', gc.master]
    coregargs = ['null', '-d', os.path.abspath(gc.procdir), '-y', gc.batchflag, '-p', gc.polygonfile,'-f', gc.polygonfile.strip('.xy'), '-z', gc.ziplistfile, '-m', gc.master]
    unwargs = ['null', '-d', os.path.abspath(gc.procdir),'-f', gc.polygonfile.strip('.xy')]
    
    L01mk.main(argv=mkargs)
    
    with cd(gc.procdir):
       demcall = 'LiCSAR_01_mk_crop_extDEM DEM/dem_crop'
       os.system(demcall)
    
    L02co.main(argv=coregargs)
    
    make_baselines(gc.procdir, mdate, gc.ziplistfile)
    
    make_interferograms_list(gc.procdir, mdate, bperpfile)
    
    unwrap_all(unwargs)
    
    
if __name__ == "__main__":
    sys.exit(main())
