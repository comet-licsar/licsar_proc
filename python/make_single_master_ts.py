#!/usr/bin/env python
""" Makes a single master timeseries from LiCSAR results
"""

import sys
import getopt
import os
import shutil
import subprocess as subp
import h5py as h5
import numpy as np
import datetime as dt

import pdb

class Usage(Exception):
    def __init__(self, msg):
        self.msg = msg
class cd:
    """Context manager for changing the current working directory"""
    def __init__(self, newPath):
        self.newPath = os.path.expanduser(newPath)

    def __enter__(self):
        self.savedPath = os.getcwd()
        os.chdir(self.newPath)

    def __exit__(self, etype, value, traceback):
        os.chdir(self.savedPath)




def main(argv=None):
    if argv == None:
        argv = sys.argv
    
    try:
        try:
            opts, args = getopt.getopt(argv[1:], "hd:s:", ["help"])
        except getopt.error, msg:
            raise Usage(msg)
        for o, a in opts:
            if o == '-h' or o == '--help':
                print __doc__
                return 0
            elif o == '-d':
                datadir = a
        if not datadir:
            raise Usage('No data directory given, -d option is not optional!')
        if not os.path.exists(datadir):
            raise Usage('Data directory {0} does not seem to exist?'.format(datadir))
    except Usage, err:
        print >>sys.stderr, "\nWoops, something went wrong:"
        print >>sys.stderr, "  "+str(err.msg)
        print >>sys.stderr, "\nFor help, use -h or --help.\n"
        return 2
    for f in os.listdir(os.path.join(datadir,'geo')):
        if f[-3:] == '.lt':
            masterdate = np.int32(f.split('.')[0])
    slavelist = []
    for l in os.listdir(os.path.join(datadir,'RSLC')):
        if l != str(masterdate) and l[0] == '2':
            slavelist.append(dt.datetime(int(l[:4]),int(l[4:6]),int(l[6:])))
    slavelist.sort()
    masterdate = dt.datetime(int(str(masterdate)[:4]),int(str(masterdate)[4:6]),int(str(masterdate)[6:]))
    
    mastermlipar = os.path.join(datadir,'SLC',masterdate.strftime('%Y%m%d'),'{md}.slc.mli.par'.format(md=masterdate.strftime('%Y%m%d')))
    mliwidth = np.int(subp.check_output(['grep','range_samples:',mastermlipar]).split(':')[1].strip())
    masterpar = os.path.join(datadir,'RSLC',masterdate.strftime('%Y%m%d'),'{md}.rslc.par'.format(md=masterdate.strftime('%Y%m%d')))
    ifgdir = os.path.join(datadir,'IFG')
    for s in slavelist:
        slavepar = os.path.join(datadir,'RSLC',s.strftime('%Y%m%d'),'{sd}.rslc.par'.format(sd=s.strftime('%Y%m%d')))
        
        ifgname = '{0}_{1}'.format(masterdate.strftime('%Y%m%d'),s.strftime('%Y%m%d'))
        create_offset(masterpar,slavepar,os.path.join(ifgdir,ifgname+'.off'),1,1)
        make_ifg(datadir,masterdate.strftime('%Y%m%d'),s.strftime('%Y%m%d'),mliwidth)

def create_offset(slc1par,slc2par,offfile,rl,al):
    exec_str = 'create_offset {slc1} {slc2} {of} 1 {rl} {al} 0'.format(slc1=slc1par,
                                                                       slc2=slc2par,
                                                                       of=offfile,
                                                                       rl=rl,
                                                                       al=al)
    os.system(exec_str)


def make_ifg(datadir,masterdate,slavedate,mliwidth):
    slcdir = os.path.join(datadir,'SLC')
    rslcdir = os.path.join(datadir,'RSLC')
    geodir = os.path.join(datadir,'geo')
    ifgdir = os.path.join(datadir,'IFG')
    if not os.path.exists(ifgdir):
        os.mkdir(ifgdir)
    exe_str = 'phase_sim_orb {sd}/{md}/{md}.slc.par {rd}/{sld}/{sld}.rslc.par '.format(sd=slcdir,
                                                                                      rd=rslcdir,
                                                                                      md=masterdate,
                                                                                      sld=slavedate)
    exe_str += '{ifd}/{md}_{sld}.off {gd}/{md}.hgt {ifd}/{md}_{sld}.sim_unw '.format(rd=rslcdir,
                                                                                   md=masterdate,
                                                                                   sld=slavedate,
                                                                                   gd=geodir,
                                                                                   ifd=ifgdir)
    exe_str += '{sd}/{md}/{md}.slc.par - - 1 1'.format(sd=slcdir,
                                                       md=masterdate)
    os.system(exe_str)
    exe_str = 'SLC_diff_intf {sd}/{md}/{md}.slc {rd}/{sld}/{sld}.rslc '.format(sd=slcdir,
                                                                               md=masterdate,
                                                                               rd=rslcdir,
                                                                               sld=slavedate)
    exe_str += '{sd}/{md}/{md}.slc.par {rd}/{sld}/{sld}.rslc.par '.format(sd=slcdir,
                                                                          md=masterdate,
                                                                          rd=rslcdir,
                                                                          sld=slavedate)

    offfile = '{ifd}/{md}_{sld}.off'.format(ifd=ifgdir,
                                            md=masterdate,
                                            sld=slavedate)
        
    exe_str += '{of} {ifd}/{md}_{sld}.sim_unw '.format(of=offfile,
                                                       md=masterdate,
                                                       sld=slavedate,
                                                       ifd=ifgdir)
    exe_str += '{ifd}/{md}_{sld}.diff 1 1 0 0 0.2 1 1'.format(ifd=ifgdir,
                                                              md=masterdate,
                                                              sld=slavedate)

    os.system(exe_str)

    exe_str = 'rasmph_pwr {ifd}/{md}_{sld}.diff {sd}/{md}/{md}.slc.mli {mw} 1 1 0 10 10'.format(ifd=ifgdir,
                                                                                md=masterdate,
                                                                                sld=slavedate,
                                                                                sd=slcdir,
                                                                                mw=mliwidth)
    os.system(exe_str)

if __name__ == "__main__":
    sys.exit(main())
