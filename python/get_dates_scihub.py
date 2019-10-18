#!/usr/bin/env python

# this will generate list of DATES (are you sure you want files?)
# over a given frame with start and end dates provided
# e.g. get_dates_scihub.py 079D_05607_131313 20141001 20160609

import s1data
import sys

frame=sys.argv[1]
startdate=sys.argv[2]
enddate=sys.argv[3]

files = s1data.get_images_for_frame(frame, str(startdate), str(enddate))

dates = []
for a in files:
    datum = a.split('_')[5].split('T')[0]
    dates.append(datum)

dateset = list(set(dates))
dateset.sort()

for date in dateset:
    print(date)
