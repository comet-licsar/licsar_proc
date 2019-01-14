#!/usr/bin/env python
"""
MySQL database query wrappers for LiCSAR processor
"""
## History ##
# 2017/11/21: DM added batch support via global variable 'glob_batch' to skip any db interaction (in the do query function)

#############

import os
import sys
import itertools
import matplotlib.path as mpltPath
import numpy as np
#import pymysql						##
#from ConfigParser impodbquerrt SafeConfigParser
import datetime as dt
#import pdb

# Local imports
import global_config as gc
from LiCSAR_lib.LiCSAR_misc import Usage
#from dbfunctions import Conn_db			##
from . import slczip

def loadziplist(ziplistfile):
    # Load full paths to the zip files into a list
    #print("Ziplist file: %s" % ziplistfile)
    with open(ziplistfile) as zf:
            ziplist = zf.read().splitlines()
    
    if len(ziplist) <= 0:
        #print "Empty ziplist file! Abort."
        raise Usage('Empty zip list file: %s' % (ziplistfile))
    else:
        return ziplist
    
def loadpolygon(polygonfile):
    # Load full paths to the zip files into a list
    #print("Polygon file: %s" % polygonfile)
    poly = np.loadtxt(polygonfile, delimiter=' ')
    if ( poly.shape[0] <= 3 ) or ( poly.shape[1] != 2 ) :
        #print "Empty (<=3 points) polygon file! Abort."
        raise Usage('Empty (<=3 points) polygon file (or bad shape, ncols!=2): %s' % (polygonfile))
    else:
        if ( poly[-1,0] != poly[0,0] ) or ( poly[-1,1] != poly[0,1] ):
            poly2 = np.append(poly, poly[0,:])
        else:
            poly2 = poly
        return poly2 # flip left-right to adopt (lat,long) convention

def get_zip_dates(ziplist):
    # Retrieve the a list of unique dates from the filenames
    dates = []
    for line in ziplist:
        date = line.split('/')[-1].split('_')[5].split('T')[0]
        if dt.date(int(date[:4]),int(date[4:6]),int(date[6:8])) in dates:
            pass
        else:
            dates.append(dt.date(int(date[:4]),int(date[4:6]),int(date[6:8])))
    #if ( len(dates) == 1 ):
    #    raise Usage('Sorry, we can work only with more than one dates now')
    if ( len(dates) < 1 ):
        raise Usage('Bad format ziplist : %s' % (ziplistfile))
    dates.sort()
    return dates

def get_burst_centers(ziplist):
    # Retrieve a unique list of burst centers and anxids
    nswath = len(slczip.SLCzip(ziplist[0]).get_burst_ANXID())
    centers = []
    anxids = []
    for line in ziplist: # loop over the zip files
        slc = slczip.SLCzip(line)       # use the zipslc class
        swaths = slc.get_swathid()
        #print swaths
        anxid_list = slc.get_burst_ANXID()
        #print anxid_list
        _,coords_list = slc.get_burstGeo_coords()  ### conventioon here is (lat,long)
        
        for swath in range(len(swaths)): # loop over the subswaths (3 for IW)
            sws = "_"+swaths[swath]
            #print sws
            for i in range(len(anxid_list[swath])): # loop over the bursts within the subswath
                anxid = anxid_list[swath][i]
                if ( str(anxid)+sws in anxids ) or ( str(anxid+1)+sws in anxids ) or ( str(anxid-1)+sws in anxids ):
                    pass
                else:
                    #print str(anxid)+sws
                    anxids.append(str(anxid)+sws)
                    centers.append(coords_list[swath][i]) # retrieve the burst centre for the same burst
    
    # Reverse the coordinate order to lon,lat:
    centers2 = []
    for coord in centers:
        centers2.append((coord[1],coord[0]))
    return anxids, centers2
    
#def get_burst_centers_all(ziplist):
    #anx, cxy = get_burst_centers(ziplist)
    #result = swaths2all(anx, cxy)
    #return result[0], (result[1], result[2])

def np2tuplearray(data):
    result=[]
    for i in range(data.shape[0]):
        result.append(tuple((data[i,0],data[i,1])))
    return result


#def swaths2all(anxids, xys ):
    #outids = []
    #outcx = []
    #outcy = []
    #for i in xrange(len(anxids)):
        #for j in xrange(len(anxids[i])):
            #if anxids[i][j]:
                #outids.append(anxids[i][j])
                #outcx.append(xys[i][j][0])
                #outcy.append(xys[i][j][1])
    #return outids, outcx, outcy

class dbquery:
    burstidlist = []
    def __init__(self,ziplistfile, polygonfile):
        self.ziplistfile = ziplistfile
        self.polygonfile = polygonfile
        self.ziplist = loadziplist(ziplistfile)
        self.datelist = get_zip_dates(self.ziplist)
        self.polygon = loadpolygon(polygonfile)
        self.allburstanxids, self.allburstcoords = get_burst_centers(self.ziplist)
        self.anxids = self.__get_bursts_in_polygon()
        
        ## Debugging tools:
        #self.bursts, self.anxids = self.get_bursts_in_polygon()
        
    def do_query(query, commit=False):
        print("Something has attempted to run a database query in batch mode... This is probably a bug, but will carry on regardless")
        return True

    def connection_established(self):
        return True

    def check_frame(self, frame):
        # checks if frame exists in database
        return True
    
    def get_dates(self):
        return get_zip_dates(self.ziplist)

    def get_bursts_in_frame(self, frame):
        # takes frame, returns list with burstid, centre_lon and 
        # centre_lat of all bursts in frame
        return self.get_bursts_in_polygon()
    
    def get_bursts_in_polygon(self):
        return [(s,[]) for s in self.anxids]
        
    def __get_bursts_in_polygon(self):            
        # Based on https://stackoverflow.com/questions/36399381/whats-the-fastest-way-of-checking-if-a-point-is-inside-a-polygon-in-python
        #allburstanxids, allburstcoords = get_burst_centers(self.ziplist)
        path = mpltPath.Path(np2tuplearray(self.polygon))
        inside2 = path.contains_points(np.asarray(self.allburstcoords))
        centers = list(itertools.compress(self.allburstcoords,inside2))
        anxids = list(itertools.compress(self.allburstanxids,inside2))
        return anxids
    
    def get_frame_bursts_on_date(self,frame,date):
        # takes frame and datetime.date object and bursts id and center coords
        # for all all bursts within the frame that were acquired on that date
        dateform = date.strftime('%Y%m%d')
        dateziplist = [s for s in self.ziplist if dateform in s]
        anxids_iw, centers = get_burst_centers(dateziplist)
        anxids = [ s.split('_')[0] for s in anxids_iw ]
        result = []
        # Reassign burst ids to nearby ones stored in self.allburstanxids
        for i in range(len(anxids)):
            anxid = int(anxids[i])
            sws = '_'+anxids_iw[i].split('_')[1]
            if str(anxid)+sws in self.anxids:
                result.append((str(anxid)+sws,[]))
            elif str(anxid+1)+sws in self.anxids:
                result.append((str(anxid+1)+sws,[]))
            elif str(anxid-1)+sws in self.anxids:
                result.append((str(anxid-1)+sws,[]))
        return result

    def get_frame_files_date(self,frame,date):
        # takes frame and one datetime.date object and returns
        # polygon name, file name and file path for all files 
        # in frame on the given date
        dateform = date.strftime('%Y%m%d')
        dateziplist = [s for s in self.ziplist if dateform in s]  
        ziplist = []
        for zf in dateziplist:
            ziplist.append((self.polygonfile, ''.join(os.path.basename(zf).split('.')[0:-1]), zf))
        return ziplist

    
    def get_burst_no(self,frame,date):
        # takes frame and datetime.date object and returns burst numbers
        # return burst_id, file and number of burst in file
        dateform = date.strftime('%Y%m%d')
        dateziplist = [s for s in self.ziplist if dateform in s]
        t,_ = get_burst_centers([dateziplist[0]])
        nswath=len(t)
        result=[]
        for anxid_iw in self.anxids:
            anxid = int(anxid_iw.split('_')[0])
            sws = '_'+anxid_iw.split('_')[1]
            for zf in dateziplist:
                zfanxids,_ = get_burst_centers([zf])
                
                if (str(anxid)+sws in zfanxids):
                    found = str(anxid)+sws
                    foundx = str(anxid)+sws
                elif (str(anxid+1)+sws in zfanxids):
                    found = str(anxid)+sws
                    foundx = str(anxid+1)+sws
                elif (str(anxid-1)+sws in zfanxids):
                    found = str(anxid)+sws
                    foundx = str(anxid-1)+sws
                else:
                    found=[]
                    continue
                
                if found:
                    swathlist = [s for s in zfanxids if sws in s]
                    n = swathlist.index(foundx)
                    rtpl = (found,os.path.basename(zf),n)
                    result.append(rtpl)
                    break
        return result
    

    
    


    ##############################################################
    # Functions from automated version...
    
    def get_frame_files_period(frame,t1,t2):
        # takes frame and two datetime.date objects and returns list returns
        # polygon name, aquisition date, file name and file path for all files 
        # in frame in the given time period
        sql_q = "select distinct polygs.polyid_name, date(files.acq_date), " \
            "files.name, files.abs_path from files " \
            "inner join files2bursts on files.fid=files2bursts.fid " \
            "inner join polygs2bursts on files2bursts.bid=polygs2bursts.bid " \
            "inner join polygs on polygs2bursts.polyid=polygs.polyid " \
            "where polygs.polyid_name='{0}' " \
            "and date(files.acq_date) between '{1}' and '{2}' "\
            "and files.pol='VV'"\
            "order by files.acq_date;".format(frame,t1,t2)
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

