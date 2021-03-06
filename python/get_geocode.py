#!/usr/bin/env python
""" Get latitude and longitude for master geometry
"""
import getopt
import sys
import os
import subprocess
import numpy as np
from LiCSAR_lib.LiCSAR_misc import *

class Usage(Exception):
    def __init__(self, msg):
        self.msg = msg

def main(argv=None):
    datadir = ''
    masterdate = ''
    if argv == None:
        argv = sys.argv
    
    try:
        try:
            opts, args = getopt.getopt(argv[1:], "hd:m:e:", ["help"])
        except getopt.error as msg:
            raise Usage(msg)
        for o, a in opts:
            if o == '-h' or o == '--help':
                print(__doc__)
                return 0
            elif o == '-d':
                datadir = a
            elif o == '-m':
                masterdate = a
    except Usage as err:
        print("\nWoops, something went wrong:", file=sys.stderr)
        print("  "+str(err.msg), file=sys.stderr)
        print("\nFor help, use -h or --help.\n", file=sys.stderr)
        return 2

    geocode_mli(datadir,masterdate)


def geocode_mli(datadir,masterdate):
    geodir = os.path.join(datadir,'geo')
    
    dempar = os.path.join(geodir,'EQA.dem_par')
    res = grep1('corner_lat',dempar)
    demlat = np.float32(res.split(':')[1].strip().split(' ')[0])
    res = grep1('corner_lon',dempar)
    demlon = np.float32(res.split(':')[1].strip().split(' ')[0])
    res = grep1('post_lat',dempar) 
    latstep = np.float32(res.split(':')[1].strip().split(' ')[0])
    res = grep1('post_lon',dempar) 
    lonstep = np.float32(res.split(':')[1].strip().split(' ')[0])
    res = grep1('width',dempar)
    dem_width = np.int32(res.split(':')[1].strip())
    res = grep1('nlines',dempar) 
    dem_length = np.int32(res.split(':')[1].strip())
    mlipar = os.path.join(datadir,'SLC',masterdate,masterdate+'.slc.mli.par')
    res = grep1('range_samples',mlipar)
    mliwidth = np.int32(res.split(':')[1].strip())
    res = grep1('azimuth_lines',mlipar)
    mlilength = np.int32(res.split(':')[1].strip())    
    
    lat = np.arange(demlat,demlat+dem_length*latstep,latstep)
    lat = lat[:dem_length]
    
    lon = np.arange(demlon,demlon+dem_width*lonstep,lonstep)
    lon = lon[:dem_width]

    LON,LAT = np.meshgrid(lon,lat)
    
    np.float32(LON).byteswap().tofile(geodir+'/lon_dem')
    np.float32(LAT).byteswap().tofile(geodir+'/lat_dem')
    exe_str = 'geocode {0}/{1}.lt_fine {0}/lon_dem {2} {0}/lon_mli {3} {4}'.format(geodir,
                                                                                   masterdate,
                                                                                   dem_width,
                                                                                   mliwidth,
                                                                                   mlilength)
    os.system(exe_str)
    exe_str = 'geocode {0}/{1}.lt_fine {0}/lat_dem {2} {0}/lat_mli {3} {4}'.format(geodir,
                                                                                   masterdate,
                                                                                   dem_width,
                                                                                   mliwidth,
                                                                                   mlilength)
    os.system(exe_str)
    return mliwidth, mlilength

#def grep(arg,file):
#    res = subprocess.check_output(['grep',arg,file])
#    return res

if __name__ == "__main__":
    sys.exit(main())
