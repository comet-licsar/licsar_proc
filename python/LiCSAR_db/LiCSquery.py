#!/usr/bin/env python
"""
MySQL database query wrappers for LiCSAR processor
"""


import os
import sys
import itertools
import pymysql
from configparser import ConfigParser
import datetime as dt
#import pdb
from numbers import Number
from shapely.geometry import Polygon
import pandas as pd

# Local imports
import global_config as gc

#get tunnel or not
parser = ConfigParser()
parser.read(gc.configfile)
#parser.get('sqlinfo','use_tunnel')
use_tunnel = bool(int(parser.get('sqlinfo','use_tunnel')))
if use_tunnel:
    from dbfunctions import Conn_tunnel_db as Conn_db
    print('warning - ssh tunnel will be used - close it when finished')
    conn, tunnel = Conn_db()
else:
    from dbfunctions import Conn_db
    #conn = Conn_db()

#if conn == 'MYSQL ERROR':
#    print('No database connection could be established')

def delete_frame_only(frame):
    sql = "delete from polygs where polygs.polyid_name='{}';".format(frame)
    print(sql)
    #res = do_query(sql, True)
    return True
    #return res

def rename_frame(frameold, framenew):
    sql = "update polygs set polyid_name='{0}' where polyid_name='{1}';".format(framenew, frameold)
    #print(sql)
    res = do_query(sql, True)
    #return True
    return res


def rename_burst(bold, bnew):
    print('WARNING, the operation of burst rename may cause inconsistencies (not tested deeply)')
    bid_old = get_bid_frombidtanx(bold)
    bid_new = get_bid_frombidtanx(bnew)
    #
    #change the burst in files2bursts
    sql = "update files2bursts set bid={0} where bid={1};".format(bid_new, bid_old)
    res = do_query(sql, True)
    #change the burst in frame definitions
    sql = "update polygs2bursts set bid={0} where bid={1};".format(bid_new, bid_old)
    res = do_query(sql, True)
    #remove bid_old from bursts
    sql = "delete from bursts where bid={0};".format(bid_old)
    res = do_query(sql, True)
    return res


def do_pd_query(query):
    if use_tunnel:
        global conn
    else:
        conn = Conn_db()
    #if type(conn) == tuple:
    #    conn = conn[0]
    df = pd.read_sql_query(query, conn)
    return df


def do_query(query, commit=False):
    # execute MySQL query and return result
    try:
        if use_tunnel:
            global conn
        else:
            conn = Conn_db()
        if conn == 'MYSQL ERROR':
            print('No database connection could be established to perform the following query: \n%s' % query)
            return 'MYSQL ERROR'
        #if type(conn) == tuple:
        #    conn = conn[0]
        #to reconnect potentially lost connection..
        rc = conn.ping(reconnect=True)
        with conn.cursor() as c:
            c.execute(query, )
            res_list = c.fetchall()
            if commit:
                conn.commit()
                res_list = c.rowcount
    except pymysql.err.Error as e:
        print("\nUnexpected MySQL error {0}: {1}".format(e[0], e[1]))
        return []
    return res_list


def close_db_and_tunnel():
    if use_tunnel:
        #print('debug - killing connection')
        try:
            conn.kill(conn.thread_id())
        except:
            print('')
        #print('debug - closing connection')
        try:
            conn.close()
        except:
            print('MySQL connection perhaps already closed?')
        #print('debug - closing tunnel')
        if tunnel.is_active:
            tunnel.close()
            print('ssh tunnel closed')
            return True
        else:
            #print('there is no ssh tunnel established here')
            return False
    else:
        try:
            conn.close()
        except:
            print('error closing db connection')
        return True


def connection_established():
    sql_q = "SELECT VERSION();"
    res = do_query(sql_q)
    if res == 'MYSQL ERROR':
        print('Could not establish database connection')
        return False
    return True


def check_frame(frame):
    # checks if frame exists in database
    sql_q = "select distinct polyid_name from polygs " \
        "where polyid_name='{0}'".format(frame)
    return do_query(sql_q)

def geom_from_polygs2geom(frame):
    polyid = get_frame_polyid(frame)[0][0]
    sql_q = "select AsText(geom) from polygs2gis where polyid={0}".format(polyid)
    return do_query(sql_q)[0][0]  #.decode('UTF-8')  --- this is needed in case of bit older version of MySQL

def get_polygon_from_bidtanx(bidtanx):
    sql = "select corner1_lat, corner1_lon, corner2_lat, corner2_lon, corner3_lat, corner3_lon," \
    " corner4_lat, corner4_lon from bursts where bid_tanx = '{0}';".format(bidtanx)
    coords = do_query(sql)[0]
    #the coordinates are 'random' so I do a loop to get valid polygon
    lat_set = [(0, 2, 4, 6),(0, 2, 6, 4),(0, 4, 6, 2),(0, 4, 2, 6),(0, 6, 4, 2),(0, 6, 2, 4)]
    lon_set = [(1, 3, 5, 7),(1, 3, 7, 5),(1, 5, 7, 3),(1, 5, 3, 7),(1, 7, 5, 3),(1, 7, 3, 5)]
    for i in range(len(lat_set)):
        lat_point_list = []
        lon_point_list = []
        for x in lat_set[i]:
            try:
                lat_point_list.append(coords[x])
            except:
                print('some error in coords, maybe mysql connection')
                return None
        for y in lon_set[i]:
            lon_point_list.append(coords[y])
        polygon_geom = Polygon(zip(lon_point_list, lat_point_list))
        if polygon_geom.is_valid:
            break
    return polygon_geom

def get_polygon_from_frame(frame):
    sql_q = "select corner1_lon, corner2_lon, corner3_lon, corner4_lon, " \
        "corner5_lon, corner6_lon, corner7_lon, corner8_lon, " \
        "corner9_lon, corner10_lon, corner11_lon, corner12_lon " \
        "from polygs where polyid_name='{0}';".format(frame)
    lons = do_query(sql_q)[0]
    sql_q = "select corner1_lat, corner2_lat, corner3_lat, corner4_lat, " \
        "corner5_lat, corner6_lat, corner7_lat, corner8_lat, " \
        "corner9_lat, corner10_lat, corner11_lat, corner12_lat " \
        "from polygs where polyid_name='{0}';".format(frame)
    lats = do_query(sql_q)[0]
    lons2 = []
    lats2 = []
    for lon in lons:
        if isinstance(lon, Number):
            lons2.append(lon)
    for lat in lats:
        if isinstance(lat, Number):
            lats2.append(lat)
    lons2.append(lons[0])
    lats2.append(lats[0])
    polygon = Polygon(zip(lons2,lats2))
    return polygon


def get_bursts_in_file(filename):
    # S1 file name can end on .zip, .SAFE or be just identifier..
    filename = filename.split('.')[0]
    sql_q = "select distinct bursts.bid_tanx from bursts " \
            "inner join files2bursts fb on fb.bid=bursts.bid " \
            "inner join files on fb.fid=files.fid "\
            "where files.name='{0}';".format(filename)
    return do_query(sql_q)


def get_bursts_in_frame(frame):
    # takes frame, returns list with burstid, centre_lon and 
    # centre_lat of all bursts in frame
    frametest = check_frame(frame)
    if not frametest:
        print('\nWarning!\nFrame {0} not found in database'.format(frame))
        return []
    else:
        sql_q = "select distinct bursts.bid_tanx, bursts.centre_lon, "\
            "bursts.centre_lat from bursts " \
            "inner join polygs2bursts on polygs2bursts.bid=bursts.bid " \
            "inner join polygs on polygs2bursts.polyid=polygs.polyid "\
            "where polygs.polyid_name='{0}';".format(frame)
        return do_query(sql_q)

def get_bidtanxs_in_frame(frame):
    # takes frame, returns list with burstid, centre_lon and 
    # centre_lat of all bursts in frame
    frametest = check_frame(frame)
    if not frametest:
        print('\nWarning!\nFrame {0} not found in database'.format(frame))
        return []
    else:
        sql_q = "select distinct bursts.bid_tanx from bursts " \
            "inner join polygs2bursts on polygs2bursts.bid=bursts.bid " \
            "inner join polygs on polygs2bursts.polyid=polygs.polyid "\
            "where polygs.polyid_name='{0}';".format(frame)
        return do_query(sql_q)


def get_bidtanxs_in_track(track = '001A', onlyFrames = True):
    # takes trackid (e.g. '001A'), returns list with burstid, centre_lon and 
    if onlyFrames:
        sql_q = "select distinct bursts.bid_tanx from bursts " \
            "inner join polygs2bursts on polygs2bursts.bid=bursts.bid " \
            "inner join polygs on polygs2bursts.polyid=polygs.polyid "\
            "where polygs.polyid_track='{0}';".format(track)
    else:
        if type(track) == str:
            track = str(int(track[:3]))
        sql_q = "select distinct bursts.bid_tanx from bursts where bid_tanxtk = {0};".format(track)
    return do_query(sql_q)


def get_frame_files_period(frame,t1,t2):
    # takes frame and two datetime.date objects and returns list returns
    # polygon name, aquisition date, file name and file path for all files 
    # in frame in the given time period
    #
    # the ordering is important here due to new versions of files
    # simply put, ESA recomputes some slcs times to times and the newer
    # (better) version is again ingested to NLA and licsinfo. We should
    # use only the newer version files.
    #
    # this cannot be: and files.abs_path not like '%metadata_only%' "\
    sql_q = "select distinct polygs.polyid_name, date(files.acq_date), " \
        "files.name, files.abs_path from files " \
        "inner join files2bursts on files.fid=files2bursts.fid " \
        "inner join polygs2bursts on files2bursts.bid=polygs2bursts.bid " \
        "inner join polygs on polygs2bursts.polyid=polygs.polyid " \
        "where polygs.polyid_name='{0}' " \
        "and date(files.acq_date) between '{1}' and '{2}' "\
        "and files.pol='VV' "\
        "order by files.acq_date asc, files.name asc, files.proc_date desc, files.date_added desc;".format(frame,t1,t2)
    return do_query(sql_q)

def get_frame_files_date(frame,date):
    # takes frame and one datetime.date object and returns
    # polygon name, file name and file path for all files 
    # in frame on the given date
    
    #in this mess, some scripts use date, and some timestamp or datetime..
    #let's convert it to date type only
    if type(date) is not type(dt.datetime.now().date()):
        date = date.date()
    
    #this is to fix for the around-midnight data:
    date2 = date + dt.timedelta(days=1)
    sql_q = "select distinct polygs.polyid_name, " \
        "files.name, files.abs_path from files " \
        "inner join files2bursts on files.fid=files2bursts.fid " \
        "inner join polygs2bursts on files2bursts.bid=polygs2bursts.bid " \
        "inner join polygs on polygs2bursts.polyid=polygs.polyid " \
        "where polygs.polyid_name='{0}' " \
        "and (date(files.acq_date)='{1}' or date(files.acq_date)='{2}')" \
        "and pol='VV'"\
        "order by files.acq_date ASC, files.date_added DESC;".format(frame,date,date2)
    return do_query(sql_q)

def get_frame_polyid(frame):
    sql_q = "select distinct polyid from polygs where polyid_name='{0}';".format(frame)
    return do_query(sql_q)

def get_framename_from_fid(fid):
    sql_q = "select distinct polyid_name from polygs where polyid={0};".format(str(fid))
    return do_query(sql_q)

def get_bid_frombidtanx(bidtanx):
    sql_q = "select distinct bid from bursts where bid_tanx='{0}';".format(bidtanx)
    return do_query(sql_q)[0][0]


def get_ipf(filename):
    # filename should be e.g. S1A_IW_SLC__1SSV_20141222T210739_20141222T210809_003837_00496E_C84D
    sql_q = "select distinct proc_vers from files where name='{0}';".format(filename)
    return do_query(sql_q)[0][0]

def get_burst_no(frame,date):
    # takes frame and datetime.date object and returns burst numbers
    # return burst_id, file and number of burst in file
    sql_q = "select distinct bursts.bid_tanx, " \
        "files.name, files2bursts.burst_no from files " \
        "inner join files2bursts on files.fid=files2bursts.fid " \
        "inner join polygs2bursts on files2bursts.bid=polygs2bursts.bid " \
        "inner join polygs on polygs2bursts.polyid=polygs.polyid " \
        "inner join bursts on polygs2bursts.bid=bursts.bid "\
        "where polygs.polyid_name='{0}' " \
        "and date(files.acq_date)='{1}' "\
        "and pol='VV'"\
        "order by files.acq_date;".format(frame,date)
    return do_query(sql_q)


def get_frame_bursts_on_date(frame,date):
    # takes frame and datetime.date object and bursts id and center coords
    # for all all bursts within the frame that were acquired on that date
    frametest = check_frame(frame)
    if not frametest:
        print('Frame {0} not found in database'.format(frame))
        return []
    else:
        sql_q = "select distinct bursts.bid_tanx, bursts.centre_lon, "\
            "bursts.centre_lat from bursts " \
            "inner join polygs2bursts on polygs2bursts.bid=bursts.bid " \
            "inner join files2bursts on files2bursts.bid=bursts.bid "\
            "inner join polygs on polygs2bursts.polyid=polygs.polyid "\
            "inner join files on files2bursts.fid=files.fid "\
            "where polygs.polyid_name='{0}'"\
            "and files.pol='VV'"\
            "and date(files.acq_date)='{1}';".format(frame,date)
        return do_query(sql_q)


def get_bursts_in_polygon_old(lon1,lon2,lat1,lat2,relorb = None, swath = None):
    #swath can be provided, e.g. as 'S6'
    #relorb can be provided, e.g. as 75
    #however this is very simplified function - if center is outside of the lat lon area, it would not find anything....
    # - but what to do if CEDA's mySQL is so historic it doesn't have geotables?
    #The GIS-enabled postgreSQL db that A.McD. was working so hard on is not used -- for ..various reasons.. but it perhaps should
    sql_q = "select distinct bursts.bid_tanx, bursts.centre_lon, bursts.centre_lat, files.rel_orb, files.swath from bursts " \
        "inner join files2bursts on files2bursts.bid=bursts.bid "\
        "inner join files on files2bursts.fid=files.fid "\
        "where bursts.centre_lon >= '{0}' and bursts.centre_lon <= '{1}' ".format(lon1,lon2)
    if swath:
        sql_q += "and files.swath='{0}' ".format(swath)
    if relorb:
        sql_q += "and files.rel_orb={0} ".format(relorb)
    sql_q += "and bursts.centre_lat >= '{0}' and bursts.centre_lat <= '{1}';".format(lat1,lat2)
    return do_query(sql_q)

def get_bursts_in_polygon(minlon,maxlon,minlat,maxlat,relorb = None, swath = None):
    if swath:
        if 'IW' not in swath:
            swath = 'IW'+str(swath)
    sql_q = "select distinct b.bid_tanx from bursts b " \
            "inner join files2bursts on files2bursts.bid=b.bid " \
            "inner join files on files2bursts.fid=files.fid " \
            "where greatest ( " \
            "b.corner1_lon, b.corner2_lon, b.corner3_lon, b.corner4_lon) >= {0} and least(" \
            "b.corner1_lon, b.corner2_lon, b.corner3_lon, b.corner4_lon) <= {1} ".format(minlon,maxlon)
    if swath:
        sql_q += "and files.swath='{0}' ".format(swath)
    if relorb:
        sql_q += "and files.rel_orb={0} ".format(relorb)
    sql_q += "and greatest( " \
             "b.corner1_lat, b.corner2_lat, b.corner3_lat, b.corner4_lat) >= {0} and least(" \
             "b.corner1_lat, b.corner2_lat, b.corner3_lat, b.corner4_lat) <= {1};".format(minlat,maxlat)
    return do_query(sql_q)

def get_frames_in_lonlat(lon,lat):
    radius = 0.01
    minlon = lon-radius
    minlat = lat-radius
    maxlon = lon+radius
    maxlat = lat+radius
    frames = get_frames_in_polygon(minlon,maxlon,minlat,maxlat)
    return frames

def get_frames_in_polygon(minlon,maxlon,minlat,maxlat):
    sql_q = "select distinct p.polyid_name from polygs p inner join polygs2bursts pb on p.polyid=pb.polyid " \
            "inner join bursts b on pb.bid=b.bid where greatest( " \
            "b.corner1_lon, b.corner2_lon, b.corner3_lon, b.corner4_lon) >= {0} and least(" \
            "b.corner1_lon, b.corner2_lon, b.corner3_lon, b.corner4_lon) <= {1} ".format(minlon,maxlon)
    sql_q += "and greatest( " \
             "b.corner1_lat, b.corner2_lat, b.corner3_lat, b.corner4_lat) >= {0} and least(" \
             "b.corner1_lat, b.corner2_lat, b.corner3_lat, b.corner4_lat) <= {1};".format(minlat,maxlat)
    return do_query(sql_q)

def get_frames_in_orbit(relorb, orbdir = None):
    relorb = str(relorb)
    while len(relorb) < 3:
        relorb = '0'+relorb
    if orbdir:
        #e.g. 124D
        relorb = relorb+orbdir
    sql_q = "select distinct polyid_name from polygs where polyid_name LIKE '{}%';".format(relorb)
    return do_query(sql_q)

def get_files_from_burst(burstid):
    sql_q = "select distinct bursts.bid_tanx, files.abs_path from bursts " \
        "inner join files2bursts on files2bursts.bid=bursts.bid "\
        "inner join files on files2bursts.fid=files.fid "\
        "where bursts.bid_tanx = '{0}';".format(burstid)
    return do_query(sql_q)


def get_filenames_from_burst(burstid):
    sql_q = "select distinct files.name from bursts " \
        "inner join files2bursts on files2bursts.bid=bursts.bid "\
        "inner join files on files2bursts.fid=files.fid "\
        "where bursts.bid_tanx = '{0}';".format(burstid)
    return do_query(sql_q)


def get_orbit_from_bidtanx(bidtanx):
    #e.g. bidtanx = '127_IW1_20509'
    sql_q = "select f.orb_dir from files f inner join files2bursts fb " \
        "on f.fid=fb.fid inner join bursts b on fb.bid=b.bid " \
        "where b.bid_tanx='{0}' limit 1;".format(bidtanx)
    return do_query(sql_q)[0][0]

def get_polygon(polyid_nm):
    sql_q = "SELECT corner1_lon, corner1_lat, corner2_lon, corner2_lat, " \
        "corner3_lon, corner3_lat, corner4_lon, corner4_lat, " \
        "corner5_lon, corner5_lat, corner6_lon, corner6_lat, " \
        "corner7_lon, corner7_lat, corner8_lon, corner8_lat, " \
        "corner9_lon, corner9_lat, corner10_lon, corner10_lat, " \
        "corner11_lon, corner11_lat, corner12_lon, corner12_lat " \
        "FROM polygs where polyid_name = '%s';" % (polyid_nm)
    return do_query(sql_q)

def get_time_of_file(fileID):
    #takes file ID like e.g. S1B_IW_SLC__1SDV_20190808T040601_20190808T040628_017489_020E3F_C558
    #and gets the full time of its acquisition, not only date..
    #you may get the file ID e.g. using
    #date=datetime.strptime('2019-08-08','%Y-%m-%d')
    #filelist = get_frame_files_date(frame, date.date())
    #fileID = filelist[0][1]
    sql_q = "SELECT acq_date from files where name ='{0}';".format(fileID)
    try:
        out = do_query(sql_q)[0][0]
    except:
        out = None
    return out

def get_wkt_boundaries(frameName):
    polygon = get_polygon(frameName)[0]
    frame_poly_cor = []
    for framepoly in polygon:
        if framepoly != None:
            frame_poly_cor.append(framepoly)
        frame_poly = frame_poly_cor
        frame_poly_zip = list(zip(frame_poly[::2], frame_poly[1::2]))
    #first is lon, second is lat
    minlon = minlat = 200
    maxlon = maxlat = -200
    for lon, lat in frame_poly_zip:
        if lon < minlon: minlon = round(lon,3)
        if lon > maxlon: maxlon = round(lon,3)
        if lat < minlat: minlat = round(lat,3)
        if lat > maxlat: maxlat = round(lat,3)
    polyText = str(minlon)+' '+str(minlat)+','+ \
               str(minlon)+' '+str(maxlat)+','+ \
               str(maxlon)+' '+str(maxlat)+','+ \
               str(maxlon)+' '+str(minlat)+','+ \
               str(minlon)+' '+str(minlat)
    wkt = 'POLYGON(({0}))'.format(polyText)
    return wkt


def is_in_polygs2geom(frameid):
    polyid = get_frame_polyid(frameid)
    if polyid:
        polyid = polyid[0][0]
    else:
        print('some error - the frame ID does not exist?')
        return None
    is_in_geom_sql = "select count(*) from polygs2gis where polyid = {0};".format(str(polyid))
    res = do_query(is_in_geom_sql)
    is_in_geom = res[0][0]
    return is_in_geom

def set_job_started(job_id):
    sql_q = "UPDATE jobs "\
        "SET licsar_version = '%s' , " \
        "  time_started = IF( "\
        "  time_started IS NULL, " \
        "  NOW(), " \
        "  time_started "\
        "), " \
        "job_status = 2 "\
        "WHERE job_id = %d;" % (gc.config['VERSION'], job_id)

    # perform query, get result (should be blank), and then commit the transaction
    res = do_query(sql_q, True)

    return


def set_job_request_started_master(job_id, acq_date):
    sql_q = "UPDATE job_requests "\
        "SET jr_status = 19 "\
        "WHERE job_id = %d " \
        "AND acq_date = DATE(%s) ;" % (job_id, acq_date)

    # perform query, get result (should be blank), and then commit the transaction
    res = do_query(sql_q, True)

    return


def set_job_request_started_standard(job_id, acq_date):
    sql_q = "UPDATE job_requests "\
        "SET jr_status = 59 "\
        "WHERE job_id = %d " \
        "AND acq_date = DATE(%s) ;" % (job_id, acq_date)

    # perform query, get result (should be blank), and then commit the transaction
    res = do_query(sql_q, True)

    return


def set_job_finished(job_id, ec):
    sql_q = "UPDATE jobs "\
        "SET time_finished = NOW(), "\
        " job_status = CASE " \
        "  WHEN %d = 0 THEN 3 " \
        "  WHEN %d > 0 AND %d < ( "\
        "    SELECT COUNT(jr_id) FROM job_requests WHERE job_id = %d "\
        "    ) "\
        "    THEN 4 " \
        "  WHEN %d > 0 AND %d = ( " \
        "    SELECT COUNT(jr_id) FROM job_requests WHERE job_id = %d " \
        "    ) " \
        "    THEN 5 " \
        "  ELSE 9 "\
        " END "\
        "WHERE job_id = %d ;" % (ec, ec, ec, job_id, ec, ec, job_id, job_id)

    # perform query, get result (should be blank), and then commit the transaction
    res = do_query(sql_q, True)

    return

def set_error_for_unclean_job_finishes(job_id):
    # If a job fails for an unknown reason and doesn't put an error code in the job_requests table, tidy it up with an
    # appropriate error code to indicate it died and the error was not caught properly.
    sql_q = "UPDATE job_requests " \
        "SET jr_status = 999 " \
        "WHERE job_id = %d " \
        "AND jr_status = 59;" % job_id

    # perform query, get result (should be blank), and then commit the transaction
    res = do_query(sql_q, True)

    return


def store_frame_geometry(frameid, wkt):
    polyid = get_frame_polyid(frameid)
    if polyid:
        polyid = polyid[0][0]
    else:
        print('some error - the frame ID does not exist?')
        return None
    is_in_geom = is_in_polygs2geom(frameid)
    if is_in_geom > 0:
        sql_q = "UPDATE polygs2gis set geom=GeomFromText('{0}') where polyid={1}".format(wkt,str(polyid))
    else:
        sql_q = "INSERT INTO polygs2gis VALUES ({0}, GeomFromText('{1}'));".format(str(polyid), wkt)
    # perform query, get result (should be blank), and then commit the transaction
    res = do_query(sql_q, True)
    return res

def sqlout2list(insql):
    #in case we get only string, we assume it was an error
    if type(insql)=='str':
        return None
    out = []
    for a in insql:
        out.append(a[0])
    return out

# DEPRECATED - delete after testing its replacement
# def set_job_finished_clean(job_id):
#     sql_q = "UPDATE jobs "\
#         "SET time_finished = NOW(), "\
#         " job_status = 3 "\
#         "WHERE job_id = %d;" % (job_id)
#
#     # perform query, get result (should be blank), and then commit the transaction
#     res = do_query(sql_q)
#     conn.commit()
#
#     return


def set_job_request_finished_clean(job_id, acq_date):
    sql_q = "UPDATE job_requests "\
        "SET jr_status = 0 "\
        "WHERE job_id = %d " \
        "AND acq_date = DATE(%s);" % (job_id, acq_date)

    # perform query, get result (should be blank), and then commit the transaction
    res = do_query(sql_q, True)

    return


# DEPRECATED - delete after testing its replacement
# def set_job_finished_error(job_id, status=9):
#     sql_q = "UPDATE jobs "\
#         "SET time_finished = NOW(), "\
#         " job_status = %d "\
#         "WHERE job_id = %d;" % (status, job_id)
#
#     # perform query, get result (should be blank), and then commit the transaction
#     res = do_query(sql_q)
#     conn.commit()
#
#     return
    

def set_job_request_finished_master_fail(job_id, acq_date, status=101):
    sql_q = "UPDATE job_requests "\
        "SET jr_status = %d "\
        "WHERE job_id = %d " \
        "AND acq_date = DATE(%s) ;" % (status, job_id, acq_date)

    # perform query, get result (should be blank), and then commit the transaction
    res = do_query(sql_q, True)

    return


def set_job_request_finished_standard_fail(job_id, acq_date, status=501):
    sql_q = "UPDATE job_requests "\
        "SET jr_status = %d "\
        "WHERE job_id = %d " \
        "AND acq_date = DATE(%s);" % (status, job_id, acq_date)

    # perform query, get result (should be blank), and then commit the transaction
    res = do_query(sql_q, True)

    return


def set_files2jobs(job_id, file_list):
    # check data types and convert it to a list of strings as required, log errors of unexepcted data types
    # these errors should only be reported if our code is supplying variable data types.
    if type(file_list) is tuple:
        if type(list(file_list)[0]) is tuple:
            file_list = [list(i) for i in file_list]
            file_list = file_list[0]
        elif type(list(file_list)) is str:
            file_list = list(file_list)

    if type(file_list) is list:
        if type(file_list[0]) is str:
            ziplist = ','.join('"{0}"'.format (f) for f in file_list)
        else:
            print("ERROR: Expected list of strings, but got list of %s" % type(file_list[0]))
    else:
        print("ERROR: Expected list of string, but got %s" % type(file_list))


    sql_q = "INSERT INTO files2jobs "\
        "(fid, job_id) "\
        "SELECT fid, %d "\
        "FROM files "\
        "WHERE abs_path "\
        "IN (%s); " % (job_id, ziplist)

    # perform query, get result (should be blank), and then commit the transaction
    res = do_query(sql_q, True)

    return


def set_new_rslc_product(job_id, acq_date, master_rslc_id, filename, filepath, rslc_status=0):
    if master_rslc_id == -1:
        sql_q = "INSERT INTO rslc "\
            "(polyid, filename, filepath, job_id, acq_date, master_rslc_id, rslc_status) "\
            "SELECT polyid, '%s', '%s', %d, date(%s), %d, %d "\
            "FROM jobs WHERE job_id = %d;"% ( filename, filepath+'/'+filename, job_id, acq_date, master_rslc_id,
                                              rslc_status, job_id)
    else:
        sql_q = "INSERT INTO rslc " \
                "(polyid, filename, filepath, job_id, acq_date, master_rslc_id, rslc_status) " \
                "SELECT j.polyid, '%s', '%s', %d, date(%s), rm.rslc_id, %d " \
                "FROM jobs  AS j "\
                "  JOIN rslc AS rm" \
                "   ON rm.filepath='%s' AND rm.rslc_status=0 "\
                "WHERE j.job_id = %d;" % (filename, filepath + '/' + filename, job_id, acq_date, rslc_status,
                                        master_rslc_id, job_id)

    # perform query, get result (should be blank), and then commit the transaction
    res = do_query(sql_q, True)

    return


def set_new_ifg_product(job_id, rslc_path_1, rslc_path_2, filepath, ifg_status=0):
    filename = filepath.split('/')[-1]
    sql_q = "INSERT INTO ifg "\
            "    (polyid, rslc_id_1, acq_date_1, rslc_id_2, acq_date_2, date_gap, filename, filepath, job_id, ifg_status) "\
            "SELECT "\
            "    j.polyid, "\
            "    r1.rslc_id, "\
            "    r1.acq_date, "\
            "    r2.rslc_id, "\
            "    r2.acq_date, "\
            "    DATEDIFF(r2.acq_date, r1.acq_date), "\
            "    '%s', "\
            "    '%s', "\
            "    %d, "\
            "    %d "\
            "FROM jobs as j "\
            "    JOIN rslc as r1 ON "\
            "        r1.filepath='%s' "\
            "    AND "\
            "        r1.rslc_status=0 "\
            "    JOIN rslc as r2 ON "\
            "        r2.filepath='%s' "\
            "    AND "\
            "        r2.rslc_status=0 "\
            "WHERE j.job_id=%d; " % (filename, filepath, job_id, ifg_status, rslc_path_1, rslc_path_2, job_id)

    # perform query, get result (should be blank), and then commit the transaction
    res = do_query(sql_q, True)

    return

def update_ifg_product_unwrapped(job_id, filename, status=1):
    sql_q = "UPDATE ifg SET unwrapped=%d WHERE job_id=%d AND filename='%s';" % (status, job_id, filename)
    # perform query, get result (should be blank), and then commit the transaction
    res = do_query(sql_q, True)
    return


def update(table='eq2frame', col='coifg_status', value='1', condition='fid=1'):
    sql_q = "UPDATE {0} SET {1}={2} WHERE {3};".format(table, col, value, condition)
    # perform query, get result (should be blank), and then commit the transaction
    res = do_query(sql_q, True)
    return


def set_new_coherence_product(job_id, rslc_path_1, rslc_path_2, filepath, coh_status=0):
    filename = filepath.split('/')[-1]
    sql_q = "INSERT INTO coherence " \
            "    (polyid, rslc_id_1, acq_date_1, rslc_id_2, acq_date_2, date_gap, filename, filepath, job_id, coh_status) " \
            "SELECT " \
            "    j.polyid, " \
            "    r1.rslc_id, " \
            "    r1.acq_date, " \
            "    r2.rslc_id, " \
            "    r2.acq_date, " \
            "    DATEDIFF(r2.acq_date, r1.acq_date), " \
            "    '%s', " \
            "    '%s', " \
            "    %d, " \
            "    %d " \
            "FROM jobs as j " \
            "    JOIN rslc as r1 ON " \
            "        r1.filepath='%s' " \
            "    AND " \
            "        r1.rslc_status=0 " \
            "    JOIN rslc as r2 ON " \
            "        r2.filepath='%s' " \
            "    AND " \
            "        r2.rslc_status=0 " \
            "WHERE j.job_id=%d; " % (filename, filepath, job_id, coh_status, rslc_path_1, rslc_path_2, job_id)

    # perform query, get result (should be blank), and then commit the transaction
    res = do_query(sql_q, True)

    return


def get_eqid(eventid):
    sql_q = "select eqid from eq where USGS_ID='{0}';".format(eventid)
    res = do_query(sql_q)
    if res:
        res = res[0][0]
    else:
        res = None
    return res


def get_usgsid(eqid):
    sql_q = "select USGS_ID from eq where eqid={0};".format(eqid)
    res = do_query(sql_q)
    if res:
        res = res[0][0]
    else:
        res = None
    return res


def insert_new_eq(event, active = True):
    if active:
        stract = '1'
    else:
        stract = '0'
    sql_q = "INSERT INTO eq " \
            "    (USGS_ID, magnitude, location, depth, time, lat, lon, active) " \
            "VALUES ('{0}',{1},'{2}',{3},'{4}',{5}, {6}, {7});".format(event.id, event.magnitude, event.location.replace("'"," "),
            event.depth, event.time.strftime('%Y-%m-%d %H:%M:%S'), round(event.latitude,2), round(event.longitude,2), stract)
    # perform query, get result (should be blank), and then commit the transaction
    res = do_query(sql_q, True)
    test = get_eqid(event.id)
    if not test:
        print('some error happened - eq was not saved to database')
        return False
    else:
        return test


def insert_new_eq2frame(eqid, fid, post_acq = False, active = True):
    if active:
        stract = '1'
    else:
        stract = '0'
    if post_acq:
        next_acq = post_acq + dt.timedelta(days=6)
        last_acq = post_acq + dt.timedelta(days=24)
        sql_q = "INSERT INTO eq2frame " \
            "    (eqid, fid, frame_status, post_acq, coifg_status, next_acq, last_acq, postifg_status, active) " \
            "VALUES ({0},{1},1,'{2}',1,'{3}','{4}', 1, {5});".format(eqid, fid, post_acq.strftime('%Y-%m-%d %H:%M:%S'), 
            next_acq.strftime('%Y-%m-%d %H:%M:%S'), last_acq.strftime('%Y-%m-%d %H:%M:%S'), stract)
    else:
        sql_q = "INSERT INTO eq2frame " \
        "    (eqid, fid)" \
        "VALUES ({0},{1});".format(eqid, fid)
    # perform query, get result (should be blank), and then commit the transaction
    res = do_query(sql_q, True)
    return res
