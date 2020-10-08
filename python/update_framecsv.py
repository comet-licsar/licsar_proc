#!/usr/bin/env python
import os, sys
import getopt

import framecare as fc

class Usage(Exception):
    """Usage context manager"""
    def __init__(self, msg):
        self.msg = msg

def main(argv=None):
    if argv == None:
        argv = sys.argv
    try:
        try:
            opts, args = getopt.getopt(argv[1:], "f:", ["version", "help"])
        except getopt.error as msg:
            raise Usage(msg)
        for o, a in opts:
            if o == '-f':
                frame = a

    except Usage as err:
        print("\nERROR:", file=sys.stderr)
        print("  "+str(err.msg), file=sys.stderr)
        print("\nFor help, use -h or --help.\n", file=sys.stderr)
        return 2

    #currentdir = '/gws/nopw/j04/nceo_geohazards_vol1/projects/LiCS/proc/current'
    #pubdir = '/gws/nopw/j04/nceo_geohazards_vol1/public/LiCSAR_products/'
    currentdir = os.environ['LiCSAR_procdir']
    pubdir = os.environ['LiCSAR_public']
    track = str(int(frame[0:3]))
    trackdir = os.path.join(currentdir,track)
    framedir = os.path.join(trackdir,frame)
    pubframedir = os.path.join(pubdir,track,frame,'metadata')
    if not os.path.exists(pubframedir):
        print('error - this frame does not exist in LiCSAR_public, cancelling')
        return False
    gpan = fc.frame2geopandas(frame)
    fc.export_frames_to_licsar_csv(gpan)
    #if rc > 0:
    #    print('updated in framecsv')
if __name__ == "__main__":
    sys.exit(main())

