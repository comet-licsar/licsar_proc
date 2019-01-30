#!/usr/bin/env python
import os
import pdb
import sys
import subprocess as subp
from shutil import copyfile

LiCSARpath = '/home/users/karzzeh/Software/LiCSAR'
sys.path.append(LiCSARpath)
sys.path.append(LiCSARpath+'/bin')
sys.path.append(LiCSARpath+'/lib')
sys.path.append(LiCSARpath+'/LiCSdb')
sys.path.append(LiCSARpath+'/python')

currentdir = '/gws/nopw/j04/nceo_geohazards_vol1/projects/LiCS/proc/current'
tracklist = [t for t in os.listdir(currentdir) if os.path.isdir(os.path.join(currentdir,t))]
pubdir = '/gws/nopw/j04/nceo_geohazards_vol1/public/LiCSAR_products/'

for track in tracklist:
    trackdir = os.path.join(currentdir,track)
    framelist = [f for f in os.listdir(trackdir) if os.path.isdir(os.path.join(trackdir,f))]
    for frame in framelist:
        if os.path.exists(os.path.join(trackdir,frame,'IFG')):
            ret = os.popen('bsub -q short-serial -W 0:05 -o {0}.log -e {0}.err "./submit_lookangles.py -t {1} -f {0}"'.format(frame,track))
