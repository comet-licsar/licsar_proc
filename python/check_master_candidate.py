#!/usr/bin/env python

# this will check if given epoch can be ok candidate for frame
# e.g. check_master_candidate.py 079D_05607_131313 20190501

import s1data
import sys
import LiCSquery as lq
import datetime as dt
from LiCSAR_lib.mk_imag_lib import check_master_bursts

framename=sys.argv[1]
m=sys.argv[2]
masterdate = dt.date(int(m[:4]),int(m[4:6]),int(m[6:8]))


print('rechecking S1 data for this date')
todown = s1data.check_and_import_to_licsinfo(framename,masterdate - dt.timedelta(days=1), masterdate + dt.timedelta(days=1), reingest = False)
filelist = lq.get_frame_files_period(framename,masterdate-dt.timedelta(days=1),masterdate+dt.timedelta(days=1))
midnighterror = False
for f in filelist:
    if f[1] > masterdate: midnighterror = True
burstlist = lq.get_bursts_in_frame(framename)
rc = check_master_bursts(framename,burstlist,masterdate,[masterdate],lq, midnighterror)

if rc != 0:
    print('ok')
else:
    print('missing bursts')
