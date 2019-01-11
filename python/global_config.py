import os
import sys
import subprocess as subp

config = {
    'VERSION': 'V2.1',
    'ENV': 'production', #[production|development|mirror]
    'DEST': 'CEMS' #[CEMS|leeds]
}

rglks = 20
azlks = 4
outres = 0.001
batchflag=False

configfile = os.environ["LiCSARconfig"]
