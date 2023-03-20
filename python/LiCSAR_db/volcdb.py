#!/usr/bin/env python
"""
MySQL database query wrappers for volcanoes
"""

from LiCSquery import *

'''
describe volc_frame_clips;
+--------------+---------------------+------+-----+---------+-------+
| Field        | Type                | Null | Key | Default | Extra |
+--------------+---------------------+------+-----+---------+-------+
| vid          | int(11)             | NO   | PRI | NULL    |       |
| volc_id      | int(8) unsigned     | NO   | MUL | NULL    |       |
| resolution_m | tinyint(3) unsigned | YES  |     | NULL    |       |
| diameter_km  | tinyint(3) unsigned | YES  |     | NULL    |       |
| polyid       | int(11)             | NO   | MUL | NULL    |       |
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

describe volc_frame_subtitutes;
+----------------+---------+------+-----+---------+-------+
| Field          | Type    | Null | Key | Default | Extra |
+----------------+---------+------+-----+---------+-------+
| vid            | int(11) | NO   | MUL | NULL    |       |
| vid_substitute | int(11) | NO   | MUL | NULL    |       |
+----------------+---------+------+-----+---------+-------+

'''

def get_volclip_info(vid=None): #,extended=True):
    """ This will load info about volcanic frame clip.
    if vid==None: it will return table for all existing clip definitions.
    """
    if vid:
        cond = "where vid={}".format(str(vid))
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

