import os
import sys
import subprocess as subp

config = {
    'VERSION': 'V2.1',
    'ENV': 'production', #[production|development|mirror]
    'DEST': 'CEMS' #[CEMS|leeds]
}

#rglks and azlks values for TOPS mode S-1 images (dest. resolution: ~100x100 m)
rglks = 20
azlks = 4
#for stripmap images (destination resolution: ~30x30 m)
rglks_stripmap = 11
azlks_stripmap = 7
#output resolution in degrees, used e.g. for DEM generation
outres = 0.001
batchflag=False

configfile = os.environ["LiCSARconfig"]

#this trick here makes variables available
#just save them to local_config.py file in your processing folder..
#cdir = os.getcwd()
sys.path.append('')
if os.path.exists('local_config.py'):
    from local_config import *
    #import local_config
