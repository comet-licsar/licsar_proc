#!/usr/bin/env python
"""
MySQL database query wrappers for LiCSAR processor
"""


import os
import sys
import itertools
import pymysql
from configparser import SafeConfigParser
import datetime as dt
import pdb

# Local imports
import global_config as gc
from dbfunctions import Conn_db

def do_query(query, commit=False):
    # execute MySQL query and return result
    try:
        conn = Conn_db()
        if conn == 'MYSQL ERROR':
            print('No database connection could be established to perform the following query: \n%s' % query)
            return 'MYSQL ERROR'

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
        "order by files.acq_date asc, files.name asc, files.proc_date desc;".format(frame,t1,t2)
    return do_query(sql_q)

def get_frame_files_date(frame,date):
    # takes frame and one datetime.date object and returns
    # polygon name, file name and file path for all files 
    # in frame on the given date
    sql_q = "select distinct polygs.polyid_name, " \
        "files.name, files.abs_path from files " \
        "inner join files2bursts on files.fid=files2bursts.fid " \
        "inner join polygs2bursts on files2bursts.bid=polygs2bursts.bid " \
        "inner join polygs on polygs2bursts.polyid=polygs.polyid " \
        "where polygs.polyid_name='{0}' " \
        "and date(files.acq_date)='{1}' "\
        "and pol='VV'"\
        "order by files.acq_date;".format(frame,date)
    return do_query(sql_q)

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

def get_bursts_in_polygon(lon1,lon2,lat1,lat2):
    sql_q = "select distinct bursts.bid_tanx, bursts.centre_lon, bursts.centre_lat, files.rel_orb, files.swath from bursts " \
        "inner join files2bursts on files2bursts.bid=bursts.bid "\
        "inner join files on files2bursts.fid=files.fid "\
        "where bursts.centre_lon >= '{0}' and bursts.centre_lon <= '{1}' ".format(lon1,lon2)
    sql_q += "and bursts.centre_lat >= '{0}' and bursts.centre_lat <= '{1}';".format(lat1,lat2)
    return do_query(sql_q)

def get_files_from_burst(burstid):
    sql_q = "select distinct bursts.bid_tanx, files.abs_path from bursts " \
        "inner join files2bursts on files2bursts.bid=bursts.bid "\
        "inner join files on files2bursts.fid=files.fid "\
        "where bursts.bid_tanx = '{0}';".format(burstid)
    return do_query(sql_q)

def get_polygon(polyid_nm):
    sql_q = "SELECT corner1_lon, corner1_lat, corner2_lon, corner2_lat, " \
        "corner3_lon, corner3_lat, corner4_lon, corner4_lat, " \
        "corner5_lon, corner5_lat, corner6_lon, corner6_lat, " \
        "corner7_lon, corner7_lat, corner8_lon, corner8_lat, " \
        "corner9_lon, corner9_lat, corner10_lon, corner10_lat, " \
        "corner11_lon, corner11_lat, corner12_lon, corner12_lat " \
        "FROM polygs where polyid_name = '%s';" % (polyid_nm)
    return do_query(sql_q)
    
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

