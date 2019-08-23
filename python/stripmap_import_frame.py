#!/usr/bin/env python
from LiCSAR_lib.LiCSAR_misc import get_centre_from_latlon, get_colat10
import LiCSAR_db.LiCSquery as lq
from batchDBLib import get_polyid
import sys, os, glob

dir = sys.argv[1]
#dir = '/gws/nopw/j04/nceo_geohazards_vol1/projects/LiCS/proc/current/stripmap/075A_SM_FOGO_S6'
if not os.path.exists(dir):
    print('wrong dir - please correct')
    print('run this script with a (full) path to the stripmap frame folder')
    exit()

xyfile = os.path.join(dir,'frame.xy')
frame = dir.split('/')[-1]
trackstr = frame.split('_')[0]
track = int(trackstr[0:-1])
swath = frame.split('_')[-1]
masterdate = os.listdir(os.path.join(dir,'SLC'))
if len(masterdate) != 1:
    print('error - you should have nothing else in SLC directory than the master dir!')
    exit()
masterdate=masterdate[0]
masterfile=glob.glob(os.path.join(dir,'SLC',masterdate,'*zip'))
if len(masterfile) != 1:
    print('error - you should have exactly one slc zip in the master SLC directory!')
    exit()
masterfile=masterfile[0].split('/')[-1]

polyid = get_polyid(frame)
if polyid:
    print('This frame has been already imported!')
    exit()

#reading through coordinates
#we assume 4 coordinates only!
coords = 4
latlon = []
with open(xyfile) as f:
    for i in range(coords):
        line = f.readline()
        lat = float(line.split(' ')[0])
        lon = float(line.split(' ')[1])
        latlon.append((lat, lon))

#get center coords in order to get colat10
centre = get_centre_from_latlon(latlon)
colat10 = get_colat10(centre[0])

#get master fid and related bursts
sql_q = "select fid from files where name='{0}' and pol='VV';".format(masterfile.split('.')[0])
sqlout = lq.do_query(sql_q)
if len(sqlout) < 1:
    print('some error happened - the master image is not ingested to database')
    exit()
if len(sqlout) > 1:
    print('some error happened - the master image appears multiple time in database')
    print('check with Milan..')
    exit()
fid = sqlout[0][0]

#get burst related to the master
sql_q = "select bid from files2bursts where fid={0};".format(fid)
sqlout = lq.do_query(sql_q)
if len(sqlout) != 1:
    print('some error happened - there should be exactly one burst ID related to the frame')
    print('check with Milan..')
    exit()
bid = sqlout[0][0]

#insert frame to licsinfo database
sql_q = "INSERT INTO polygs " \
        "(polyid_name, polyid_track, polyid_colat10, corner1_lat, corner1_lon, corner2_lat, corner2_lon, " \
        "corner3_lat, corner3_lon, corner4_lat, corner4_lon) " \
        "values ('{0}', '{1}', {2}, {3}, {4}, {5}, {6}, {7}, {8}, {9}, " \
        "{10});".format(frame, trackstr, colat10, latlon[0][0], latlon[0][1], latlon[1][0], latlon[1][1], \
        latlon[2][0], latlon[2][1], latlon[3][0], latlon[3][1])
sqlout = lq.do_query(sql_q, True)

if sqlout != 1:
    print('error in inserting frame def. to database')
    exit()

polyid = get_polyid(frame)

#connect frame ID with the burst ID
sql_q = "INSERT INTO polygs2bursts VALUES ({0}, {1});".format(polyid, bid)
sqlout = lq.do_query(sql_q, True)

if sqlout == 1:
    print('success, frame imported as polyid={0} and connected to burst {1}'.format(polyid, bid))
