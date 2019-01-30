#!/usr/bin/env python
import os
import glob
import pdb
import sys
import subprocess as subp
from shutil import copyfile, copy
import getopt
import numpy as np


LiCSARpath = '/home/users/karzzeh/Software/LiCSAR'
sys.path.append(os.getcwd())
sys.path.append(LiCSARpath)
sys.path.append(LiCSARpath+'/bin')
sys.path.append(LiCSARpath+'/lib')
sys.path.append(LiCSARpath+'/LiCSdb')
sys.path.append(LiCSARpath+'/python')

from LiCSAR_02_coreg import get_mli_size, get_dem_size
from gamma_functions_karzzeh import look_vector, geocode

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
    pubframedir = os.path.join(pubdir,track,frame,'metadata')
    if not os.path.exists(os.path.join(framedir,'log')):
        os.mkdir(os.path.join(framedir,'log'))
    origmasterdate = os.listdir(os.path.join(framedir,'SLC'))[0]
    geodir = os.path.join(framedir,'geo')
    slcdir = os.path.join(framedir,'SLC',origmasterdate)
    slcpar = os.path.join(slcdir,origmasterdate+'.slc.par')
    [width, length] = get_mli_size(os.path.join(slcdir,origmasterdate+'.slc.mli.par'))
    offpar = '-'
    dempar = os.path.join(framedir,'geo','EQA.dem_par')
    [demwidth, demlength] = get_dem_size(os.path.join(dempar))
    dem = os.path.join(framedir,'geo','EQA.dem')
    theta = os.path.join(geodir,'theta')
    phi = os.path.join(geodir,'phi')
    lutfile = os.path.join(geodir,origmasterdate+'.lt_fine')
    logfilename = os.path.join(framedir,'log','look_vector-{0}.log'.format(origmasterdate))
    if look_vector(slcpar,offpar,dempar,dem,theta,phi,logfilename):
        logfilename = os.path.join(framedir,'log','geocode-theta-{0}.log'.format(origmasterdate))
        geocode(lutfile,theta,demwidth,theta+'.rc',width,length,'0','0',logfilename)
        logfilename = os.path.join(framedir,'log','geocode-phi-{0}.log'.format(origmasterdate))
        geocode(lutfile,phi,demwidth,phi+'.rc',width,length,'0','0',logfilename)
        thetarc = np.fromfile(theta+'.rc',dtype=np.float32).byteswap().reshape((int(length),int(width)))
        nanix = thetarc == 0
        thetarc[nanix] = np.nan
        phirc = np.fromfile(phi+'.rc',dtype=np.float32).byteswap().reshape((int(length),int(width)))
        phirc[nanix] = np.nan
        U = np.sin(thetarc)
        E = np.cos(phirc)*np.cos(thetarc)
        N = np.sin(phirc)*np.cos(thetarc)
        #theta = np.fromfile(theta,dtype=np.float32).byteswap()
        #phi = np.fromfile(phi,dtype=np.float32).byteswap()
        #U = np.sin(theta)
        #E = np.cos(phi)*np.cos(theta)
        #N = np.sin(phi)*np.cos(theta)
        U[nanix] = 0
        E[nanix] = 0
        N[nanix] = 0
        U.byteswap().tofile(os.path.join(geodir,'U'))
        E.byteswap().tofile(os.path.join(geodir,'E'))
        N.byteswap().tofile(os.path.join(geodir,'N'))
        os.remove(theta)
        os.remove(theta+'.rc')
        os.remove(phi)
        os.remove(phi+'.rc')
        #os.system('create_geoctiff_lookangles.sh {0} {1}'.format(framedir,origmasterdate))
        os.system(os.getcwd()+'/create_geoctiff_lookangles.sh {0} {1}'.format(framedir,origmasterdate))
        for tif in glob.glob(os.path.join(pubframedir,'*.tif')):
            os.remove(tif)
        for tif in glob.glob(os.path.join(framedir,'GEOC','lookangles','*.tif')):
            copy(tif,pubframedir)

if __name__ == "__main__":
    sys.exit(main())

