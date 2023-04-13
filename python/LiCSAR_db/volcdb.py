#!/usr/bin/env python
"""
MySQL database query wrappers for volcanoes
ML 2023
"""

from LiCSquery import *
from shapely.geometry import Polygon
from dbfunctions import Conn_sqlalchemy
import geopandas as gpd
import framecare as fc
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


def get_volcanoes_in_polygon(polygon, volcs = None):
    """Will get volcanoes in the gpd (e.g., volcs=get_volc_info() ) that are within given polygon (shapely.geometry.Polygon).
    You may create the polygon from lon/lat borders, e.g. using fc.lonlat_to_poly(lon1, lon2, lat1, lat2)"""
    if type(volcs) == type(None):
        volcs = get_volc_info()
    isin = volcs.geometry.apply(lambda x: polygon.contains(x))
    return volcs[isin]


def vis_gpd(vgpd):
    """Simple example to plot volcanoes or volclips in gpd.DataFrame"""
    fc.vis_aoi(vgpd.geometry.to_list())


def find_volcano_by_name(name, volcstable = None):
    """This will try to find the volcano in volcs table extracted using get_volc_info(). If the table is not given, it will extract it.
    
    Args:
        name (str): e.g. 'Askja'
        volcstable (pd.DataFrame): must contain column 'name'
    
    Returns:
        record of found volc
    """
    if type(volcstable) == type(None):
        volcstable = get_volc_info()
    return volcstable[volcstable.name.str.contains(name)]


def is_in_volclips(volcid):
    """Checks if the volcano (volcid) has its volclip."""
    return is_in_table(volcid, 'volc_id', 'volclip2volcs')


def create_volclip_for_volcano(volcid, cliparea_geo = None):
    """Will create a volclip definition for given volcano (volcid).
    If cliparea is not set, it will auto-generate area using diameter 25 km centered on the volc lon/lat.
    
    Args:
        volcid (int): the volcano ID
        cliparea_geo (str): OPTIONAL: clip boundaries, e.g. 'lon1/lon2/lat1/lat2', but will accept also shapely.geometry.Polygon (!!)
    """
    if is_in_volclips(volcid):
        print('ERROR, this volcano has already its volclip, cancelling for now')
        print('you can delete the volclip for the volcano using: delete_volclip(vid)')
        print('where you can get vid using: get_volclip_vids(volcid)')
        return False
    #
    vpd = get_volc_info(volcid)
    if vpd.empty:
        print('no record found for this volcano ID - please check again')
        return False
    #
    if type(cliparea_geo) == type(None):
        # prepare polygon
        clon, clat = vpd.lon.values[0], vpd.lat.values[0]
        radius_km = 25/2
        radius_deg=radius_km/111
        lon1=clon-radius_deg
        lon2=clon+radius_deg
        lat1=clat-radius_deg
        lat2=clat+radius_deg
    elif type(cliparea_geo) == type('string'):
        from lics_unwrap import cliparea_geo2coords
        lon1, lon2, lat1, lat2 = cliparea_geo2coords(cliparea_geo)
    else:
        lon1,lat1,lon2,lat2=cliparea_geo.bounds
        lon1,lon2=sorted([lon1,lon2])
        lat1,lat2=sorted([lat1,lat2])
    #
    lonlats = [(lon1,lat1), (lon1,lat2), (lon2,lat2), (lon2,lat1), (lon1,lat1)]
    polygon = Polygon(lonlats)
    wkt = polygon.wkt
    #
    sql_q="SELECT MAX(vid) from volclips;"
    lastvid=do_query(sql_q)[0][0]
    vid=lastvid+1  # because of auto-increment
    #
    # adding to volclips
    sql_q = "INSERT INTO volclips (vid, geometry) VALUES ({0}, GeomFromText('{1}'));".format(str(vid), wkt)
    res = do_query(sql_q, True)
    #
    #sql_q = "select last_insert_id();"
    # link the vid and volcano:
    sql_q = "INSERT INTO volclip2volcs (vid, volc_id) VALUES ({0}, {1});".format(str(vid), str(volcid))
    res = do_query(sql_q, True)
    return vid


def add_volcano_to_volclip(volcid, vid):
    """This will add a 'volcid' volcano to a 'vid' volclip"""
    sql_q = "INSERT INTO volclip2volcs (vid, volc_id) VALUES ({0}, {1});".format(str(vid), str(volcid))
    res = do_query(sql_q, True)
    return res


def delete_volclip(vid):
    """This will delete volcano clip definition - careful, not many checks in place
    """
    sql_q = "DELETE FROM volclip2volcs where vid = {0};".format(str(vid))
    res = do_query(sql_q, True)
    sql_q = "DELETE FROM volclips where vid = {0};".format(str(vid))
    res = do_query(sql_q, True)
    return res


'''
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
'''

def init_all_subsets():
    """This will auto-init all volclips (assuming only one vid per volcano...)
    """
    volcs=get_volc_info()
    for i,volc in volcs.iterrows():
        print(volc['name'])
        frames = get_volcano_frames(volc['volc_id'])
        vid = get_volclip_vids(volcid)[0]
        if frames:
            for frame in frames:
                initialise_subset_volclip(vid, frame)


def get_volcano_frames(volcid):
    """Gets frames covering the given volcano.
    """
    volc=get_volc_info(volcid)
    #print(volc['name'])
    try:
        frames= lq.sqlout2list(lq.get_frames_in_lonlat(volc.lon,volc.lat))
    except:
        print('no frame found')
        frames = False
    return frames


def get_volclip_vids(volcid):
    """Gets all volclip vids for given volcano.
    volcid is the volcano ID"""
    sql = "select vf.vid from volclip2volcs vf inner join volcanoes v on vf.volc_id=v.volc_id where v.volc_id={};".format(str(volcid))
    a=do_query(sql)
    a=sqlout2list(a)
    return a


def get_volcano_from_vid(vid):
    sql = "select volc_id from volclip2volcs where vid={};".format(str(vid))
    a=do_query(sql)
    a=sqlout2list(a)[0]
    return a


def get_volclips_gpd(vid=None):
    """Gets volclips as geodatabase - either one if given vid, or all"""
    if vid:
        cond = " where vc.vid={}".format(str(vid))
    else:
        cond = ''
    #sql = "SELECT ST_AsBinary(geometry) as geom from volclips {0};".format(cond)
    sql = "SELECT v.volc_id,v.name,vc.vid,ST_AsBinary(vc.geometry) as geom from volclips vc inner join volclip2volcs vf on vf.vid=vc.vid inner join volcanoes v on vf.volc_id=v.volc_id {0};".format(cond)
    engine=Conn_sqlalchemy()
    volclips = gpd.GeoDataFrame.from_postgis(sql, engine, geom_col='geom' )
    return volclips


def initialise_subset_volclip(vid, frame = None, resol_m = 30):
    """This will initialise the volclip subset for given frame and volclip.
    If frame is 'None', it will do this for all relevant (fully overlapping) frames.
    """
    vclipdb = get_volclips_gpd(vid)
    #volcid = vclipdb.volc_id.values[0]
    lon1, lat1, lon2, lat2 = list(vclipdb.geom.bounds.values[0])
    if not frame:
        print('getting related frames')
        frames = fc.subset_get_frames(lon1, lon2, lat1, lat2, full_overlap=True, only_initialised=True)
        for frame in frames:
            #fc.subset_initialise_corners(frame, lon1, lon2, lat1, lat2, sid = str(volcid), is_volc = True, resol_m=resol_m)
            fc.subset_initialise_corners(frame, lon1, lon2, lat1, lat2, sid = str(vid), is_volc = True, resol_m=resol_m)
    else:
        #fc.subset_initialise_corners(frame, lon1, lon2, lat1, lat2, sid = str(volcid), is_volc = True, resol_m=resol_m)
        fc.subset_initialise_corners(frame, lon1, lon2, lat1, lat2, sid = str(vid), is_volc = True, resol_m=resol_m)


# HOW TO MERGE:
'''
# select volcanoes in the area (in Kamchatka):
lat1,lon1=56.275212713142146, 161.08682836667546
lat2,lon2=55.67535919777727, 160.1041122831717
polygon=lonlat_to_poly(lon1, lon2, lat1, lat2)
selclips=get_volcanoes_in_polygon(polygon, volclips)
vis_gpd(selclips)

# say that we will merge only last 5 clips of the selection
tomerge=selclips.tail(5)
lon1,lat1,lon2,lat2=tomerge.geom.cascaded_union.bounds
poly=lonlat_to_poly(lon1, lon2, lat1, lat2)
fc.vis_subset_frames(lon1, lon2, lat1, lat2)

# we now want to remove all vids (volclips) and instead map the volcids to the new vid.
# so let's just choose one volcid, create the volclip for it, and then attach the other ones:
vids=tomerge.vid.values
volcids=tomerge.volc_id.values
volcid1=volcids[0]
# delete the older (assuming unused!) vids:
for vid in vids:
    delete_volclip(vid)
# create new volclip with volcid1
newvid = create_volclip_for_volcano(volcid1, cliparea_geo = poly)
# attach other volcanoes to this volclip
for volcid in volcids[1:]:
    add_volcano_to_volclip(volcid, newvid)

# ok, and now finally, initialise that volclip
frames = subset_get_frames(lon1, lon2, lat1, lat2, full_overlap=True, only_initialised=True)
for frame in frames:
    fc.subset_initialise_corners(frame, lon1, lon2, lat1, lat2, sid = str(vid), is_volc = True)

# and if you want, start their processing:
for frame in frames:
    cmd = 'framebatch_update_frame.sh -P -u '+frame+' upfill'
    os.system(cmd)
'''

'''
frames=['147A_02466_191712','009D_02504_202119','045A_02494_171816']
for frame in frames:
    vdb.initialise_subset_volclip(vid, frame)

frame='111D_02490_152021'
vid = 2712


volcid = get_volcano_from_vid(vid)
    
    147A_02466_191712
    111D_02490_152021
    009D_02504_202119
    045A_02494_171816
'''
    
