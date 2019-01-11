#!/usr/bin/env python
"""
A script that does all the final tidy ups, such as write the time finished to the DB, and set any
closing status codes in the DB.
"""


import os
import sys
import getopt
import datetime as dt

dir_path = os.path.dirname(os.path.realpath(__file__))
sys.path.append(dir_path[:-4])
sys.path.append(dir_path[:-4]+'/bin')
sys.path.append(dir_path[:-4]+'/lib')
sys.path.append(dir_path[:-4]+'/LiCSdb')
sys.path.append(dir_path[:-4]+'/python')

import global_config as gc
import LiCSquery as lq
from gamma_functions import *


class Usage(Exception):
    """Usage context manager"""
    def __init__(self, msg):
        self.msg = msg

def main(argv=None):
    if argv == None:
        argv = sys.argv

    framename = []
    procdir = []
    error_count = -1
    job_id = -1

    try:
        try:
            opts, args = getopt.getopt(argv[1:], "vhf:d:j:e:", ["version", "help"])
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
            elif o == '-j':
                job_id = int(a)
            elif o == '-e':
                error_count = int(a)
        if not framename:
            raise Usage('No frame, polygon or set of burst ids given, please use either the -f')
        if not procdir:
            raise Usage('No output data directory given, -d is not optional!')
        if job_id == -1:
            print "This processing is not outputting any products to the database."
        else:
            if error_count == -1:
                raise Usage('Cannot have a job_id, and no error count (-e), please supply a valid -e int')
    except Usage, err:
        print >>sys.stderr, "\nERROR:"
        print >>sys.stderr, "  "+str(err.msg)
        print >>sys.stderr, "\nFor help, use -h or --help.\n"

        return 2

    #Check if a DB connection can be established
    if not lq.connection_established():
        print >> sys.stderr, "\nERROR:"
        print >> sys.stderr, "Could not establish a stable database connection. No processing can happen."

        return 1

    # If there is a job id, then determine number of error jobs from requests, and if it is a complete or partial fail
    if job_id != -1:
        lq.set_job_finished(job_id, error_count)

        lq.set_error_for_unclean_job_finishes(job_id)



if __name__ == "__main__":
    sys.exit(main())
