import ConfigParser
import sys
import os.path
import urllib


######Set up variables in model.cfg before using
dpath = os.path.dirname(__file__)
config = ConfigParser.RawConfigParser()
config.read('%s/model.cfg'%(dpath))


def ECMWFdload(bdate,hr,filedir,humidity='Q'):

        from .ecmwfapi import ECMWFDataServer

        url = "https://api.ecmwf.int/v1" 
        emid = config.get('ECMWF', 'email')
        key = config.get('ECMWF', 'key')
        server = ECMWFDataServer(url=url, key=key, email=emid)

        assert humidity in ('Q','R'), 'Unknown humidity field for ECMWF'
        if humidity in 'Q':
                humidparam = 133
        elif humidity in 'R':
                humidparam = 157
        flist = []

        for k in xrange(len(bdate)):
                day = bdate[k]
                fname = '%s/ERA-Int_%s_%s.grb'%(filedir,day,hr)
                print 'Downloading %d of %d: %s '%(k+1,len(bdate),fname)

                flist.append(fname)    
                indict = {'dataset'  : "interim",
                          'date'     : "%s"%(day),
                          'time'     : "%s"%(hr),
                          'step'     : "0",
                          'levtype'  : "pl",
                          'levelist' : "all",
                          'type'     : "an",
                          'grid'     : "128",
                          'param'    : "129/130/%3d"%(humidparam),
                          'target'   : "%s"%(fname)}

                print indict

                if not os.path.exists(fname):
                        server.retrieve(indict)
                        

        return flist





def ECMWF_olddload(bdate,hr,filedir,humidity='Q'):

        from ecmwf import ECMWFDataServer
        print 'ECMWF server has been updated. Use new server settings.'

        emid = config.get('ECMWF_old','email')
        key = config.get('ECMWF_old','key')
        server = ECMWFDataServer(
               'http://data-portal.ecmwf.int/data/d/dataserver/',      
               key, emid)


        assert humidity in ('Q','R'), 'Unknown humidity field for ECMWF'
        if humidity in 'Q':
                humidparam = 133
        elif humidity in 'R':
                humidparam = 157
        flist = []

        for k in xrange(len(bdate)):
                day = bdate[k]
                fname = '%s/ERA-Int_%s_%s.grb'%(filedir,day,hr)
                print 'Downloading %d of %d: %s '%(k+1,len(bdate),fname)

                flist.append(fname)     
                if not os.path.exists(fname):
                        server.retrieve({
                          'dataset'  : "interim_full_daily",
                          'date'     : "%s"%(day),
                          'time'     : "%s"%(hr),
                          'step'     : "0",
                          'levtype'  : "pl",
                          'levelist' : "all",
                          'type'     : "an",
                          'grid'     : "128",
                          'param'    : "129/130/%3d"%(humidparam),
                          'target'   : "%s"%(fname),
                        })

        return flist


def NARRdload(bdate,hr,filedir):

        flist = []      
        for k in xrange(len(bdate)):
                day = bdate[k]
                webdir = day[0:6]
                fname = 'narr-a_221_%s_%s00_000.grb'%(day,hr)
                flist.append('%s/%s'%(filedir,fname))
                weburl='http://nomads.ncdc.noaa.gov/data/narr/%s/%s/%s'%(webdir,day,fname)
                dname = '%s/%s'%(filedir,fname)
                print 'Downloading %d of %d: %s'%(k+1,len(bdate),fname)
                if not os.path.exists(dname):
                        urllib.urlretrieve(weburl,dname) #,reporthook)

        return flist

def ERAdload(bdate,hr,filedir):
        import cookielib
        import urllib2

        emid = config.get('ERA','email')
        pwd = config.get('ERA','key')

        cj = cookielib.CookieJar()
        opener = urllib2.build_opener(urllib2.HTTPCookieProcessor(cj))
        login_data = urllib.urlencode({'action': 'login', 'passwd': pwd, 'email' : emid})
        opener.open('https://rda.ucar.edu/cgi-bin/login',login_data)

        flist = []
        rlist = []
        for k in xrange(len(bdate)):
                day = bdate[k]
                webdir = day[0:6]
                fname = 'ei.oper.an.pl/%s/ei.oper.an.pl.regn128sc.%s%s'%(webdir,day,hr)
                flist.append(fname)
                url = 'http://rda.ucar.edu/data/ds627.0/%s'%(fname)
                resp = opener.open(url)
                fout = '%s/ERA_%s_%s.grb'%(filedir,day,hr)
                rlist.append(fout)
                print 'Downloading %d of %d'%(k+1,len(bdate))
                if not os.path.exists(fout):
                        resp = opener.open(url)
                        fid = open(fout,'w')
                        fid.write(resp.read())
                        fid.close()

        return rlist
                                
def MERRAdload(bdate,hr,filedir, hdf=True):
    flist = []
    for i in xrange(len(bdate)):
        date = bdate[i]
        filename = '%s/merra-%s-%s.hdf' %(filedir,date,hr)
        flist.append(filename)
        print 'Downloading %d of %d: %s'%((i+1),len(bdate),filename)
        yr = date[0:4]
        mon = date[4:6]
        hr = hr
        year = int(yr)
        url1 = 'http://goldsmr3.sci.gsfc.nasa.gov/daac-bin/OTF/HTTP_services.cgi?FILENAME=%2Fdata%2Fs4pa%2FMERRA%2FMAI6NPANA.5.2.0%2F'
        url2 = '%2F'
        url3 = '%2FMERRA200.prod.assim.inst6_3d_ana_Np.'
        url3n = '%2FMERRA300.prod.assim.inst6_3d_ana_Np.'
        url3o = '%2FMERRA100.prod.assim.inst6_3d_ana_Np.'

        if hdf:
            url4 = '.hdf&FORMAT=SERGLw&BBOX=-90%2C-180%2C90%2C180&TIME=1979-01-01T'
        else:
            url4 = '.hdf&FORMAT=TmV0Q0RGLw&BBOX=-90%2C-180%2C90%2C180&TIME=1979-01-01T'

        url5 = '%3A00%3A00%2F1979-01-01T'
        url6 = '%3A00%3A00&LABEL=MERRA200.prod.assim.inst6_3d_ana_Np.'
        url6n = '%3A00%3A00&LABEL=MERRA300.prod.assim.inst6_3d_ana_Np.'
        url6o = '%3A00%3A00&LABEL=MERRA100.prod.assim.inst6_3d_ana_Np.'
        if hdf:
            url7 = '.SUB.hdf&FLAGS=&SHORTNAME=MAI6NPANA&SERVICE=SUBSET_LATS4D&LAYERS=&VERSION=1.02&VARIABLES=ps%2Ch%2Ct%2Cqv'
        else:
            url7 = '.SUB.nc&FLAGS=&SHORTNAME=MAI6NPANA&SERVICE=SUBSET_LATS4D&LAYERS=&VERSION=1.02&VARIABLES=ps%2Ch%2Ct%2Cqv'

        if year < 1993:
            weburl = '%s%s%s%s%s%s%s%s%s%s%s%s%s' %(url1,yr,url2,mon,url3o,date,url4,hr,url5,hr,url6o,date,url7)
        elif year < 2001:
            weburl = '%s%s%s%s%s%s%s%s%s%s%s%s%s' %(url1,yr,url2,mon,url3,date,url4,hr,url5,hr,url6,date,url7)
        else:
            weburl = '%s%s%s%s%s%s%s%s%s%s%s%s%s' %(url1,yr,url2,mon,url3n,date,url4,hr,url5,hr,url6n,date,url7)
        dir = '%s' %(filename)
        if not os.path.exists(dir):
            urllib.urlretrieve(weburl,dir)

    return flist


############################################################
# Program is part of PyAPS                                 #
# Copyright 2012, by the California Institute of Technology#
# Contact: earthdef@gps.caltech.edu                        #
############################################################
