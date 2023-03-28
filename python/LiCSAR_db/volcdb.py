#!/usr/bin/env python
"""
MySQL database query wrappers for volcanoes
ML 2023
"""

from LiCSquery import *
from shapely.geometry import Polygon
from dbfunctions import Conn_sqlalchemy
import geopandas as gpd
gpd.io.file.fiona.drvsupport.supported_drivers['KML'] = 'rw'

'''
CREATE TABLE volclips (vid INT(11) UNSIGNED AUTO_INCREMENT PRIMARY KEY NOT NULL, geometry POLYGON NOT NULL);
describe volclips;
+--------------+---------------------+------+-----+---------+-------+
| Field        | Type                | Null | Key | Default | Extra |
+--------------+---------------------+------+-----+---------+-------+
| vid          | int(11) unsigned    | NO   | PRI | NULL    | auto_increment |
| geometry     | polygon             | NO   |     | NULL    |       |
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

CREATE TABLE volclip2volcs (vid INT(11) UNSIGNED NOT NULL, volc_id INT(8) UNSIGNED NOT NULL);
# the 'REFERENCES' does nothing anymore:
# CREATE TABLE volclip2volcs (vid INT(11) UNSIGNED NOT NULL REFERENCES volclips(vid), volc_id INT(8) UNSIGNED NOT NULL REFERENCES volcanoes(volc_id));
describe volclip2volcs;
+----------------+---------+------+-----+---------+-------+
| Field          | Type    | Null | Key | Default | Extra |
+----------------+---------+------+-----+---------+-------+
| vid            | int(11) | NO   |     | NULL    |       |
| volc_id        | int(11) | NO   |     | NULL    |       |
+----------------+---------+------+-----+---------+-------+

'''

"""use example:
import volcdb as vdb
res = vdb.get_volc_info()
res
"""

"""
how to export the volclips to a kml?
just:
vv=get_volclips_gpd()

""" 

def export_all_volclips_to_kml(outkml):
    """e.g. outkml='/gws/nopw/j04/nceo_geohazards_vol1/public/shared/temp/earmla/volclips.kml'"""
    volclips=get_volclips_gpd()
    volclips.to_file(outkml, driver='KML')


def get_volc_info(volcid=None):
    """"This will get info on all volcanoes from the db.
    If volcid is provided, it will show info only for that volcano.
    
    try e.g. volcid=357070
    """
    if volcid:
        cond = " where volc_id={}".format(str(volcid))
    else:
        cond = ''
    sql = "select volc_id,name,lat,lon,alt,priority, ST_AsText(geometry) as geometry from volcanoes"+cond+";"
    a = do_pd_query(sql)
    a['geometry'] = a.geometry.apply(wkt.loads)
    return a


def is_in_volclips(volcid):
    """Checks if the volcano (volcid) has its volclip."""
    return is_in_table(volcid, 'volc_id', 'volclip2volcs')


def create_volclip_for_volcano(volcid):
    """Will create a volclip definition for given volcano"""
    if is_in_volclips(volcid):
        print('ERROR, this volcano has already its volclip, cancelling for now')
        return False
    
    vpd = get_volc_info(volcid)
    if vpd.empty:
        print('no record found for this volcano ID - please check again')
        return False
    
    # prepare polygon
    clon, clat = vpd.lon.values[0], vpd.lat.values[0]
    radius_km = 25/2
    radius_deg=radius_km/111
    lon1=clon-radius_deg
    lon2=clon+radius_deg
    lat1=clat-radius_deg
    lat2=clat+radius_deg
    lonlats = [(lon1,lat1), (lon1,lat2), (lon2,lat2), (lon2,lat1), (lon1,lat1)]
    polygon = Polygon(lonlats)
    wkt = polygon.wkt
    
    sql_q="SELECT MAX(vid) from volclips;"
    lastvid=do_query(sql_q)[0][0]
    # adding to volclips
    sql_q = "INSERT INTO volclips(geometry) VALUES (GeomFromText('{0}'));".format(wkt)
    res = do_query(sql_q, True)
    vid=lastvid+1  # because of auto-increment
    #sql_q = "select last_insert_id();"
    
    # link the vid and volcano:
    sql_q = "INSERT INTO volclip2volcs (vid, volc_id) VALUES ({0}, {1});".format(str(vid), str(volcid))
    res = do_query(sql_q, True)
    return vid



def get_volclip_info(vid=None): #,extended=True):
    """NOT WORKING AS POLYID ARE DROPPED FOR NOW"""
    """ This will load info about volcanic frame clip.
    if vid==None: it will return table for all existing clip definitions.
    """
    if vid:
        cond = " where vid={}".format(str(vid))
    else:
        cond = ''
    sql = "select vf.*,v.name,v.lat,v.lon,v.alt,v.priority,p.polyid_name from volclips vf \
        inner join volcanoes v on vf.volc_id=v.volc_id inner join polygs p on vf.polyid=p.polyid"+cond+";"
    a=do_pd_query(sql)
    if a.empty:
        return False
    else:
        return a


def get_volclip_vids(volcid):
    """Gets all volclip vids for given volcano.
    volcid is the volcano ID"""
    sql = "select vf.vid from volclip2volcs vf inner join volcanoes v on vf.volc_id=v.volc_id where v.volc_id={};".format(str(volcid))
    a=do_query(sql)
    a=sqlout2list(a)
    return a


def get_volclips_gpd(vid=None):
    """Gets volclips as geodatabase - either one if given vid, or all"""
    if vid:
        cond = " where vid={}".format(str(vid))
    else:
        cond = ''
    sql = "SELECT ST_AsBinary(geometry) as geom from volclips {0};".format(cond)
    sql = "SELECT v.volc_id,v.name,vc.vid,ST_AsBinary(vc.geometry) as geom from volclips vc inner join volclip2volcs vf on vf.vid=vc.vid inner join volcanoes v on vf.volc_id=v.volc_id {0};".format(cond)
    engine=Conn_sqlalchemy()
    volclips = gpd.GeoDataFrame.from_postgis(sql, engine, geom_col='geom' )
    return volclips
