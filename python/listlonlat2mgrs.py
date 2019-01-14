#!/usr/bin/env python

import sys
import mgrs
import numpy

file=sys.argv[1]
lonlat = numpy.loadtxt(file)
nburst = len(lonlat)

for n in range (0,nburst):
    m = mgrs.MGRS()
    a = m.toMGRS(lonlat[n,1],lonlat[n,0],MGRSPrecision=1)

    print(a)
