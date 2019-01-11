#!/usr/bin/env python
import numpy as np
import scipy as sp
import os
import sys
import getopt
from ConfigParser import SafeConfigParser

import global_config as gc
import LiCSquery as lq
from LiCSAR_02_coreg import get_mli_size
from gamma_functions import *

import pdb

class Usage(Exception):
    """Usage context manager"""
    def __init__(self, msg):
        self.msg = msg

def main(argv=None):
    if argv == None:
        argv = sys.argv
    # Parameter initialisation and checking
    framename=[]
    procdir=[]

    try:
        try:
            opts, args = getopt.getopt(argv[1:], "vhf:d:", ["version", "help"])
        except getopt.error, msg:
            raise Usage(msg)
        for o, a in opts:
            if o == '-h' or o == '--help':
                print __doc__
                return 0
            elif o == '-v' or o == '--version':
                print ""
                print "Current version: %s" % gc.config['VERSION']
                print ""
                return 0
            elif o == '-f':
                framename = a
            elif o == '-d':
                procdir = a

        if not procdir:
            raise Usage('No data directory given, -d is not optional!')

    except Usage, err:
        print >>sys.stderr, "\nERROR:"
        print >>sys.stderr, "  "+str(err.msg)
        print >>sys.stderr, "\nFor help, use -h or --help.\n"
        return 2

    #Check if a DB connection can be established
#    if not lq.connection_established():
#        print >> sys.stderr, "\nERROR:"
#        print >> sys.stderr, "Could not establish a stable database connection. No processing can happen."

#        return 1

    imdates = [str(k) for k in os.listdir(os.path.join(procdir,'RSLC')) 
               if len(k) == 8 
               and k[0] == '2' 
               and os.path.isdir(os.path.join(procdir,'RSLC',k))]

    ifgdates = [str(k) for k in os.listdir(os.path.join(procdir,'IFG'))
                if len(k) == 17
                and k[0] =='2'
                and os.path.isdir(os.path.join(procdir,'IFG',k))]

    origmasterdate = [x[:8] for x in os.listdir(os.path.join(procdir,'geo')) if x[-4:] == '.hgt']
    width, length = get_mli_size(os.path.join(procdir,'SLC',origmasterdate+'.slc.mli.par'))

    uw = np.zeros((len(ifgdates),length,width),dtype=np.float16)
    A = np.zeros((len(ifgdates),len(imdates)),dtype=np.int16)
    masterlist = []
    slavelist = []

    for ix, ifgd in enumerate(ifgdates):
        ifgdir = os.path.join(procdir,'IFG',ifgd)
        uwfile = os.path.join(ifgdir,ifgd+'.unw')
        uw[ix] = np.float16(np.fromfile(uwfile,dtype=np.float32).byteswap().reshape((length,width)))
        masterdate = ifgd[:8]
        masterlist.append(masterdate)
        masterix = np.where(imdates==masterdate)[0][0]
        slavedate = ifgd[-8:]
        slavelist.append(slavedate)
        slaveis = np.where(imdates==slavedate)[0][0]
        A[ix,masterix] = -1
        A[ix,slaveix] = 1
        # Needs clever referencing here!

    A[:,0] = 0
    phvec = uw.flatten().reshape((uw.shape[0],uw.shape[1]*uw.shape[2]))
    phsm = np.zeros((len(imdates),phvec.shape[1]),dtype=np.float16)

    datedt = [dt.date(int(d[:4]),int(d[4:6]),int(d[6:])) for d in imdates]
    day = np.array([d.toordinal() for d in datedt])[:,None]
    day = day-day[0]
      
    G = np.concatenate((np.ones_like(day),day),axis=1)
    v = np.zeros_like(phsm[0,:])
    phres = np.zeros_like(phvec)*np.nan

    for ix in range(phvec.shape[1]):
        if ii%1000 == 0:
            print ix
            phvecthis= np.float32(phvec[:,ii])
        ix = np.where(np.isnan(phvecthis))[0]
        ix2 = np.where(~np.isnan(phvecthis))[0]
        if len(ix) == len(phvecthis):
            ph_sm[:,ii] = np.nan
            v[ii] = np.nan
            continue
        Athis = A.copy()
        Athis = np.delete(Athis,ix,axis=0)
        if np.linalg.matrix_rank(Athis) < Athis.shape[1]-1:
            ph_sm[:,ii] = np.nan
            v[ii] = np.nan
        else:
            phvecthis = np.delete(phvecthis,ix)
            phsmthis = np.linalg.lstsq(Athis,phvecthis)[0]
            diff = np.diff(phsmthis)
            phresthis = np.dot(Athis,phsmthis)-phvecthis 
            resix = np.where(np.abs(phresthis) > 0.8*np.pi)[0]
            # Doing crude check for unwrapping errors, rejecting suspects and re-estimating
            if len(resix) > 0:
                phvecthis = np.delete(phvecthis,resix)
                Athis = np.delete(Athis,resix,axis=0)
                if np.linalg.matrix_rank(Athis) < Athis.shape[1]-1:
                    ph_sm[:,ii] = np.nan
                    v[ii] = np.nan
                else:
                    phsmthis = np.linalg.lstsq(Athis,phvecthis)[0]
                    phresthis = np.dot(Athis,phsmthis)-phvecthis 
                    ix2 = np.delete(ix2,resix)
                    #phsmthis *= 28/2/np.pi
                    phres[ix2,ii] = phresthis
                    ph_sm[:,ii] = phsmthis
                    #dummy, vthis = np.linalg.lstsq(G,phsmthis)[0]
                    #v[ii] = vthis
            else:
                #phsmthis *= 28/2/np.pi
                phres[ix2,ii] = phresthis
                ph_sm[:,ii] = phsmthis
                #dummy, vthis = np.linalg.lstsq(G,phsmthis)[0]
                #v[ii] = vthis

        
        
    

if __name__ == "__main__":
    sys.exit(main())
