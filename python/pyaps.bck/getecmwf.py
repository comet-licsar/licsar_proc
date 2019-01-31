#!/usr/bin/python
from ecmwf import ECMWFDataServer

#####Replace with your ECMWF key and email address.
####server = ECMWFDataServer(
####       'http://data-portal.ecmwf.int/data/d/dataserver/',
####       'a9a3c89533c035e5a92c5758fbf27875',
####       'abc@def.com'
####    )


def getfiles(bdate,hr,filedir,humidity='Q')

    server = ECMWFDataServer('http://data-portal.ecmwf.int/data/d/dataserver/',
        'YOURPASSWORDATECMWF',
        'YOURLOGIN@ECMWF.INT')

    assert humidity in ('Q','R'), 'Unknown humidity field for ECMWF'
    if humidity in 'Q':
        humidparam = 133
    elif humidity in 'R':
        humidparam = 157

    for day in bdate:
        server.retrieve({
          'dataset'  : "interim_full_daily",
          'date'     : "%s"%(bdate),
          'time'     : "%s"%(hr),
          'step'     : "0",
          'levtype'  : "pl",
          'levelist' : "all",
          'type'     : "an",
          'grid'     : "128",
          'param'    : "129/130/%d"%(humidparam),
          'target'   : "%s/ERA-Int_%s_%s.grb"%(fileloc,bdate,hr),
        })


############################################################
# Program is part of PyAPS                                 #
# Copyright 2012, by the California Institute of Technology#
# Contact: earthdef@gps.caltech.edu                        #
############################################################
