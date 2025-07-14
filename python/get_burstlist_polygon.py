#!/usr/bin/env python
"""
Generate burstlist and filelist files needed for input into LiCSAR v1.0
"""

import getopt
import os
import datetime
import sys
import numpy as np
import datetime as dt
import matplotlib.path as path
from configparser import SafeConfigParser

import LiCSquery as lq 

import pdb

class Usage(Exception):
    def __init__(self, msg):
        self.msg = msg

def main(argv=None):
    if argv == None:
        argv = sys.argv
    # Parameter initialisation and checking
    polygonfile=[]
    try:
        try:
            opts, args = getopt.getopt(argv[1:], "hp:", ["help"])
        except getopt.error as msg:
            raise Usage(msg)
        for o, a in opts:
            if o == '-h' or o == '--help':
                print(__doc__)
                return 0
            elif o == '-p':
                polygonfile = a
        if not polygonfile:
            raise Usage('Polygonfile not given')
        if not os.path.exists(polygonfile):
            raise Usage('Polygon file {0} does not seem to exist'.format(polygonfile))

    except Usage as err:
        print("\nERROR:", file=sys.stderr)
        print("  "+str(err.msg), file=sys.stderr)
        print("\nFor help, use -h or --help.\n", file=sys.stderr)
        return 2
    burstidfilename = polygonfile.split('.')[0]+'_burst_ids.txt'
    filesfilename = polygonfile.split('.')[0]+'_zipfiles.list'
    with open(polygonfile,'r') as f:
        poly = np.loadtxt(f,np.float32,delimiter=' ')
        burststmp = lq.get_bursts_in_polygon(poly[:,0].min(),
                                             poly[:,0].max(),
                                             poly[:,1].min(),
                                             poly[:,1].max())
        burstids = [b[0] for b in burststmp]
        burstslon = [np.float32(b[1]) for b in burststmp] 
        burstslat = [np.float32(b[2]) for b in burststmp]
        tracks = [np.int32(b[3]) for b in burststmp]
        burstssw = [b[4] for b in burststmp]
        polypath = path.Path(poly)
        ix = polypath.contains_points(np.array((burstslon,burstslat)).T)
        tracks = [t for t,i in zip(tracks,ix) if i]
        burstids = [b for b,i in zip(burstids,ix) if i]
        burstslon = [b for b,i in zip(burstslon,ix) if i]
        burstslat = [b for b,i in zip(burstslat,ix) if i]
        burstssw = [b for b,i in zip(burstssw,ix) if i]
    print('Available tracks:\n {0}'.format(''.join([str(t)+' ' for t in set(tracks)])))


    trackchoice = np.int8(input('Please chose a track...\n'))
    if trackchoice in tracks:
        lonchoice = [b for b,t in zip(burstslon,tracks) if trackchoice == t]
        latchoice = [b for b,t in zip(burstslat,tracks) if trackchoice == t]
        idchoice = [b for b,t in zip(burstids,tracks) if trackchoice == t]
        swathchoice = [b for b,t in zip(burstssw,tracks) if trackchoice == t]
        count = []
        for sw in ['IW1','IW2','IW3']:
           count.append(len([s for s in swathchoice if s == sw]))
        with open(burstidfilename,'w') as f:
            f.write('IW1 {0} IW2 {1} IW3 {2}\n'.format(count[0],count[1],count[2]))
            for lon,lat in zip(lonchoice,latchoice):
                f.write('burstID {0} {1}\n'.format(lon,lat))
        filelist = []
        for b in idchoice:
            res = lq.get_files_from_burst(b)
            for e in res:
                filelist.append(e[1])
        filelist.sort()
        with open(filesfilename,'w') as f:
            for fl in set(filelist):
                f.write('{0}\n'.format(fl))
        pdb.set_trace()
    else:
        print('Selected track {0} not a valid option...'.format(trackchoice))

if __name__ == "__main__":
    sys.exit(main())
