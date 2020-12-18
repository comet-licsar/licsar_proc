import os
import sys
import subprocess as subp


#rglks and azlks values for TOPS mode S-1 images (dest. resolution: ~100x100 m)
rglks = 20
azlks = 4
#for stripmap images (destination resolution: ~30x30 m)
#rglks_stripmap = 11
#azlks_stripmap = 7

#output resolution in degrees, used e.g. for DEM generation
outres = 0.001
#outres = 0.0005 #this value is safe (max) to geocode ifgs as they are 20/4 multilooked

#alpha value and window used for spatial filtering by GAMMA's adf
adf_alpha = 1
adf_window = 32

#this threshold is used for coherence in unwrapping
coh_unwrap_threshold = 0.35
# ... the unwrapping is actually quite advanced! so... keeping the threshold very very low...
#coh_unwrap_threshold = 0.15
# however this causes the unwrapping step very, very long..
# so a compromise:
#coh_unwrap_threshold = 0.3


##########################################
# the lines below should not be edited
##########################################
config = {
    'VERSION': 'V2.1',
    'ENV': 'production', #[production|development|mirror]
    'DEST': 'CEMS' #[CEMS|leeds]
}

batchflag=False

configfile = os.environ["LiCSARconfig"]

#this trick here makes variables available
#just save them to local_config.py file in your processing folder..
sys.path.append('')
if os.path.exists('local_config.py'):
    from local_config import *
