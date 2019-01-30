#!/usr/bin/env python
import os
import pdb
import sys
import subprocess as subp
from shutil import copy
import getopt
import numpy as np


LiCSARpath = '/home/users/karzzeh/Software/LiCSAR'
sys.path.append(LiCSARpath)
sys.path.append(LiCSARpath+'/bin')
sys.path.append(LiCSARpath+'/lib')
sys.path.append(LiCSARpath+'/LiCSdb')
sys.path.append(LiCSARpath+'/python')

class Usage(Exception):
    """Usage context manager"""
    def __init__(self, msg):
        self.msg = msg

def main(argv=None):
    if argv == None:
        argv = sys.argv

    try:
        try:
            opts, args = getopt.getopt(argv[1:], "f:t:", ["version", "help"])
        except getopt.error as msg:
            raise Usage(msg)
        for o, a in opts:
            if o == '-f':
                frame = a
            elif o == '-t':
                track = a

    except Usage as err:
        print("\nERROR:", file=sys.stderr)
        print("  "+str(err.msg), file=sys.stderr)
        print("\nFor help, use -h or --help.\n", file=sys.stderr)
        return 2

    currentdir = '/gws/nopw/j04/nceo_geohazards_vol1/projects/LiCS/proc/current'
    pubdir = '/gws/nopw/j04/nceo_geohazards_vol1/public/LiCSAR_products/'
    trackdir = os.path.join(currentdir,track)
    framedir = os.path.join(trackdir,frame)
    framepubdir= os.path.join(pubdir,track,frame,'products')
    ifgdir = os.path.join(framedir,'IFG')
    ifglist = [x for x in os.listdir(ifgdir) if len(x) == 17]
    for ifg in ifglist:
        masterdate = ifg[:8]
        slavedate = ifg[-8:]
        basecall = 'base_orbit {0} {1} - | grep perpendicular'.format(os.path.join(framedir,'RSLC',masterdate,masterdate+'.rslc.par'),
                                                                      os.path.join(framedir,'RSLC',slavedate,slavedate+'.rslc.par'))
        try:
            bperpthis  = subp.check_output(basecall,shell=True).split(':')[1].strip()
        except:
            continue
        basefile = os.path.join(ifgdir,ifg,ifg+'.base')
        try:
            with open(basefile,'w') as f:
                f.write(bperpthis+' m\n')
            ifgpubdir = os.path.join(framepubdir,ifg)
            copy(basefile,ifgpubdir)
        except:
            continue

if __name__ == "__main__":
    sys.exit(main())

