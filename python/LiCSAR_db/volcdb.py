#!/usr/bin/env python
"""
MySQL database query wrappers for volcanoes
ML 2023
"""

from LiCSquery import *

'''
describe volc_clips;
+--------------+---------------------+------+-----+---------+-------+
| Field        | Type                | Null | Key | Default | Extra |
+--------------+---------------------+------+-----+---------+-------+
| vid          | int(11)             | NO   | PRI | NULL    |       |
| geometry     | polygon             | NO   |     |         |       |
+--------------+---------------------+------+-----+---------+-------+

describe volcanoes;
+----------+-----------------+------+-----+---------+-------+
| Field    | Type            | Null | Key | Default | Extra |
+----------+-----------------+------+-----+---------+-------+
| volc_id  | int(8) unsigned | NO   | PRI | NULL    |       |
| name     | varchar(40)     | NO   |     | NULL    |       |
| lat      | float(7,5)      | NO   |     | NULL    |       |
| lon      | float(8,5)      | NO   |     | NULL    |       |
| alt      | float(5,1)      | YES  |     | NULL    |       |
| priority | char(2)         | YES  |     | NULL    |       |
| geometry | point           | YES  |     | NULL    |       |
+----------+-----------------+------+-----+---------+-------+

describe volclip_to_volcs;
+----------------+---------+------+-----+---------+-------+
| Field          | Type    | Null | Key | Default | Extra |
+----------------+---------+------+-----+---------+-------+
| vid            | int(11) | NO   | MUL | NULL    |       |
| volc_id        | int(11) | NO   | MUL | NULL    |       |
+----------------+---------+------+-----+---------+-------+

'''

"""use example:
import volcdb as vdb
res = vdb.get_volc_info()
res
"""

def get_volc_info(volcid=None):
    """"This will get info on all volcanoes from the db.
    If volcid is provided, it will show info only for that volcano.
    """
    if volcid:
        cond = " where volc_id={}".format(str(volcid))
    else:
        cond = ''
    sql = "select volc_id,name,lat,lon,alt,priority, ST_AsText(geometry) as geometry from volcanoes"+cond+";"
    a = do_pd_query(sql)
    a['geometry'] = a.geometry.apply(wkt.loads)
    return a


def is_in_volclips(volc_id):
    return  is_in_table(volc_id, 'volc_id', 'volclip_to_volcs')


def create_volclip(volc_id, lon1, lon2, lat1, lat2):
    if is_in_volclips:
        print('ERROR, this volcano has already its volclip, cancelling for now')
        return False
    vid = #.............
    dwkt = #......... create wkt from the lons/lats polygon
    sql_q = "INSERT INTO volc_clips (vid, geometry) VALUES ({0}, GeomFromText('{1}'));".format(str(vid), dwkt)
    res = do_query(sql_q, True)
    time.sleep(0.25)
    sql_q = "INSERT INTO volclip_to_volcs (vid,volc_id) VALUES ({0}, {1});".format(str(volc_id), vid)
    res = do_query(sql_q, True)
    return vid

def get_volclip_info(vid=None): #,extended=True):
    """ This will load info about volcanic frame clip.
    if vid==None: it will return table for all existing clip definitions.
    """
    if vid:
        cond = " where vid={}".format(str(vid))
    else:
        cond = ''
    sql = "select vf.*,v.name,v.lat,v.lon,v.alt,v.priority,p.polyid_name from volc_frame_clips vf \
        inner join volcanoes v on vf.volc_id=v.volc_id inner join polygs p on vf.polyid=p.polyid"+cond+";"
    a=lq.do_pd_query(sql)
    if a.empty:
        return False
    else:
        return a


def get_volclip_vids(volcid):
    """Gets all volclip vids for given volcano.
    volcid is the volcano ID"""
    sql = "select vf.vid from volc_frame_clips vf inner join volcanoes v on vf.volc_id=v.volc_id where v.volc_id={};".format(str(volcid))
    a=lq.do_query(sql)
    a=lq.sqlout2list(a)
    return a

