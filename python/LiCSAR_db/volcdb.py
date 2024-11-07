#!/usr/bin/env python
"""
MySQL database query wrappers for volcanoes
ML 2023
"""
import pandas as pd
import glob, os
try:
    from LiCSquery import *
    from dbfunctions import Conn_sqlalchemy
    from sqlalchemy import text
except:
    print('error loading LicsInfo tools - volcdb is quite useless then')

from shapely.geometry import Polygon
import geopandas as gpd
import framecare as fc
import fiona
try:
    gpd.io.file.fiona.drvsupport.supported_drivers['KML'] = 'rw'
except:
    pass
    #print('error, export to kml not working (update geopandas)')

fiona.drvsupport.supported_drivers['KML'] = 'rw'

try:
    subvolcpath = os.path.join(os.environ['LiCSAR_procdir'], 'subsets', 'volc') #/volc/267)
except:
    print('no LiCSAR_procdir ? subvolcpath is set wrong (but leaving for now)')
    subvolcpath = 'subsets/volc'

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
| vportal_name | varchar(40)   | YES  |     | NULL    |       |
| vportal_area | varchar(20)   | YES  |     | NULL    |       |
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
    sql = "select volc_id,name,lat,lon,alt,priority,vportal_area,vportal_name,ST_AsBinary(geometry) as geom from volcanoes"+cond+";"
    engine=Conn_sqlalchemy()
    with engine.connect() as conn:
        a = gpd.GeoDataFrame.from_postgis(text(sql), conn, geom_col='geom')
    #a = do_pd_query(sql)
    #a['geometry'] = a.geometry.apply(wkt.loads)
    return a


def get_existing_rslcs_for_volclip(vid):
    volclippath = os.path.join(subvolcpath, str(vid))
    rtable = pd.DataFrame(columns=['vid','subframe', 'no_rslcs'])
    if not os.path.exists(volclippath):
        return rtable
    for subframe in os.listdir(volclippath):
        rslcs = os.listdir(os.path.join(volclippath, subframe, 'RSLC'))
        norslcs = len(rslcs)-1
        rtable.loc[len(rtable)] = [vid,subframe, norslcs]
    return rtable


def get_sourceframe_volc(vid, subframe):
    volclippath = os.path.join(subvolcpath, str(vid))
    volcsubframepath = os.path.join(volclippath, subframe)
    try:
        frame = glob.glob(volcsubframepath + '/corners_clip*')[0].split('.')[-1]
    except:
        print('error getting sourceframe for '+volcsubframepath)
        frame = ''
    return frame


def get_status_table_volc(volcid):
    vids = get_volclip_vids(volcid)
    vidtable = pd.DataFrame(columns=['volc_id','vid','subframe', 'sourceframe', 'no_rslcs'])
    for vid in vids:
        rtable = get_existing_rslcs_for_volclip(vid)
        if not rtable.empty:
            for i,r in rtable.iterrows():
                sourceframe = get_sourceframe_volc(vid, r.subframe)
                vidtable.loc[len(vidtable)] = [volcid, vid, r.subframe, sourceframe, r.no_rslcs]
    return vidtable


def get_status_details_volc(volcid = None, vid = None):
    '''Provide either volcano ID or volclip id (vid) to extract RSLCs per subframe'''
    if not vid:
        if not volcid:
            print('provide either volcid or vid')
            return
        else:
            vid = get_volclip_vids(volcid)[0]
    vtb = get_status_table_volc(volcid)
    import os
    volclip_path = os.path.join(os.environ['LiCSAR_procdir'], 'subsets', 'volc', str(vid))
    allrslcs = []
    for subframe in vtb['subframe'].values:
        rslcs = os.listdir(os.path.join(volclip_path, subframe, 'RSLC'))
        allrslcs.append(rslcs)
        #print(len(rslcs))
    return vtb.subframe.values, allrslcs


def plot_status_details_volc(volcid=None, vid=None):
    '''just copied from jupyter ntb'''
    frames, allrslcs = get_status_details_volc(volcid=volcid, vid=vid)
    import pygmt
    import numpy as np
    import pandas as pd
    import datetime as dt
    #
    today = dt.datetime.now()
    fig = pygmt.Figure()
    ys = np.ones(len(allrslcs)).cumsum()
    miny, maxy = 0, len(ys) + 1
    minx, maxx = dt.datetime(2014, 10, 1), today
    #
    fig.basemap(
        region=[minx, maxx, miny, maxy],
        projection="X12c/4c",
        frame=["xaf+ltime", "yaf+lsubframe", "WSrt+texisting hires clips"],
    )
    #
    for i in range(len(allrslcs)):
        fill = 'black'
        rslcs = allrslcs[i]

        x = pd.to_datetime(rslcs)
        y = np.ones(len(rslcs)) * ys[i]

        fig.plot(
            x=x,
            y=y,
            style="p0.05c",
            fill=fill,
            # Set the legend label,
            transparency=25,  # set transparency level for all symbols
        )
    #
    # fig.legend(transparency=30)  # set transparency level for legends
    return fig


def get_status_table_all_volcs():
    """Gets current processing status of all volcanoes (number of already existing RSLCs of their volclips).
    Will get table in the form of:
    | volc_id | vid | subframe | number of RSLC clips |
    Note that only info for already initialised volclips is returned. Also note that number of RSLCs excludes reference epoch.

    Returns:
        pandas.DataFrame
    """
    print('Getting processing status info for all volcanoes (ETA few minutes)')
    volcs = get_volc_info()
    vtable = pd.DataFrame(columns=['volc_id','vid','subframe', 'sourceframe', 'no_rslcs'])
    for volcid in volcs.volc_id:
        volcrecs = get_status_table_volc(volcid)
        if not volcrecs.empty:
            vtable = pd.concat([vtable, volcrecs])
            vtable = vtable.reset_index(drop=True)
    return vtable


def get_volcanoes_in_frame(frame):
    '''will get all volcanoes in given frame'''
    polygon = get_polygon_from_frame(frame)
    return get_volcanoes_in_polygon(polygon)


def get_volcanoes_in_polygon(polygon, volcs = None):
    """Will get volcanoes in the gpd (e.g., volcs=get_volc_info() ) that are within given polygon (shapely.geometry.Polygon).
    You may create the polygon from lon/lat borders, e.g. using fc.lonlat_to_poly(lon1, lon2, lat1, lat2)"""
    if type(volcs) == type(None):
        volcs = get_volc_info()
    isin = volcs.geom.apply(lambda x: polygon.contains(x))
    return volcs[isin]


def vis_gpd(vgpd):
    """Simple example to plot volcanoes or volclips in gpd.DataFrame"""
    fc.vis_aoi(vgpd.geom.to_list())


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


def classify_volc(volcid, toclass = ''):
    """Will classify volcano to a given class (e.g. 'C'). If empty string is provided (default), will set to None
    """
    if toclass:
        sql_q = "UPDATE volcanoes SET priority='{0}' where volc_id = {1};".format(toclass, str(volcid))
    else:
        sql_q = "UPDATE volcanoes SET priority=NULL where volc_id = {0};".format(str(volcid))
    res = do_query(sql_q, True)
    print('done')


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
    """This will delete volcano clip definition - careful, not many checks in place.
    To delete from the core system in terminal for now:
    vid=
    mv $LiCSAR_procdir/subsets/volc/$vid $LiCSAR_procdir/subsets/volc.todel/.
    for x in $LiCSAR_procdir/subsets/volc.todel/$vid/*/corner*; do
      a=`basename $x`; a=`echo $a | cut -d '.' -f2`; cdproc $a; rm subsets/$vid; cd -;
    done
    """
    sql_q = "DELETE FROM volclip2volcs where vid = {0};".format(str(vid))
    res = do_query(sql_q, True)
    sql_q = "DELETE FROM volclips where vid = {0};".format(str(vid))
    res = do_query(sql_q, True)
    return res


def replace_volclips_by_larger_area(new_cliparea_geo):
    """ This will delete db records (only) of volclips in given cliparea and link to volcanoes in this area

    Args:
        new_cliparea_geo (str):  licsbas style string, i.e. lon1/lon2/lat1/lat2

    Returns:
        gpd.DataFrame:  table for the new volclip with its volcanoes
    """
    #lon1, lon2, lat1, lat2 = -22.725, -21.6375, 63.787, 64.03855
    #new_cliparea_geo = str(lon1) + '/' + str(lon2) + '/' + str(lat1) + '/' + str(lat2)
    lon1, lon2, lat1, lat2 = np.array(new_cliparea_geo.split('/')).astype(float)
    roi = fc.lonlat_to_poly(lon1, lon2, lat1, lat2)
    volcs = get_volcanoes_in_polygon(roi)  # , volcstable)
    #
    # now delete the smaller volclips of those volcanoes
    firstrun = True
    for v in volcs.volc_id:
        # print(v)
        vids = get_volclip_vids(v)
        for vid in vids:
            # print(vid)
            delete_volclip(vid)
        if firstrun:
            newclip = create_volclip_for_volcano(v, new_cliparea_geo)
        else:
            add_volcano_to_volclip(v, newclip)
        firstrun = False
    #
    volclip = get_volclips_gpd(newclip)
    return volclip


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

def init_all_subsets(volcid = None, full_overlap=True, only_classed=True):
    """This will auto-init all volclips (assuming only one vid per volcano...):
    
    Args:
        volcid (int): if not given, will init volclips for ALL volcanoes. otherwise only for the given one
        full_overlap (bool): if False, it will skip checking for full overlap of the volclip and frames. good for the last auto-run
        only_classed (bool): if True, it will skip volcanoes that are not classified (that have 'None' as the 'priority'/class)
    """
    volcs=get_volc_info(volcid)
    if only_classed:
        volcs=volcs[~volcs.priority.isnull()]
        if volcs.empty:
            print('error, the selection is without classification. Please use classify_volc first')
            return False
    for i,volc in volcs.iterrows():
        print(volc['name'])
        frames = get_volcano_frames(int(volc['volc_id']))
        vid = get_volclip_vids(int(volc['volc_id']))[0]
        if frames:
            for frame in frames:
                try:
                    initialise_subset_volclip(vid, frame,full_overlap=full_overlap)
                except:
                    print('error with init, cleaning')
                    subsetdir = '/gws/nopw/j04/nceo_geohazards_vol1/projects/LiCS/proc/current/subsets/volc/'+str(vid)+'/'+frame[:4]
                    if os.path.exists(subsetdir):
                        os.system('rm -rf '+subsetdir)

"""
def init_volcs_in_frame(frame, full_overlap=True):
    '''will init all volcs within given frame'''
    volcs = get_volcanoes_in_frame(frame)
    for i,volc in volcs.iterrows():
        print(volc['name'])
        vid = get_volclip_vids(int(volc['volc_id']))[0]
        try:
            initialise_subset_volclip(vid, frame,full_overlap=full_overlap)
        except:
            print('error with init, cleaning')
            subsetdir = '/gws/nopw/j04/nceo_geohazards_vol1/projects/LiCS/proc/current/subsets/volc/'+str(vid)+'/'+frame[:4]
            if os.path.exists(subsetdir):
                os.system('rm -rf '+subsetdir)
"""

def get_volcano_frames(volcid):
    """Gets frames covering the given volcano.
    """
    volc=get_volc_info(volcid)
    #print(volc['name'])
    try:
        frames= sqlout2list(get_frames_in_lonlat(float(volc.lon),float(volc.lat)))
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

def get_licsbas_clipstring_volclip(vid):
    """Gets clip string of LiCSBAS from a volclip"""
    return get_licsbas_clipstring_geom(get_volclips_gpd(vid).head(1).geom)


def get_licsbas_clipstring_geom(geom):
    """Gets clip string of LiCSBAS from a record of Geoseries

    Returns:
        string (in form 'lon1/lon2/lat1/lat2' )
    """
    a = round(float(geom.bounds.minx), 5)
    b = round(float(geom.bounds.maxx), 5)
    c = round(float(geom.bounds.miny), 5)
    d = round(float(geom.bounds.maxy), 5)
    return str(a)+'/'+str(b)+'/'+str(c)+'/'+str(d)


def get_licsbas_clipstring_volcano(volcid, customradius_km = 55.55/2):
    ''' Useful to get the coords as in volcano portal (0.5 degree diameter for TOPS, 0.1 for stripmap)'''
    g = get_volc_info(volcid)['geom']
    clon, clat = g[0].x, g[0].y
    radius_deg = customradius_km / 111.111
    lon1 = round(clon - radius_deg, 5)
    lon2 = round(clon + radius_deg, 5)
    lat1 = round(clat - radius_deg, 5)
    lat2 = round(clat + radius_deg, 5)
    return str(lon1)+'/'+str(lon2)+'/'+str(lat1)+'/'+str(lat2)


def get_volclips_gpd(vid=None):
    """Gets volclips as geodatabase - either one if given vid, or all"""
    if vid:
        cond = " where vc.vid={}".format(str(vid))
    else:
        cond = ''
    #sql = "SELECT ST_AsBinary(geometry) as geom from volclips {0};".format(cond)
    sql = "SELECT v.volc_id,v.name,vc.vid,ST_AsBinary(vc.geometry) as geom from volclips vc inner join volclip2volcs vf on vf.vid=vc.vid inner join volcanoes v on vf.volc_id=v.volc_id {0};".format(cond)
    engine=Conn_sqlalchemy()
    with engine.connect() as conn:
        volclips = gpd.GeoDataFrame.from_postgis(text(sql), conn, geom_col='geom')
    #volclips = gpd.GeoDataFrame.from_postgis(sql, engine, geom_col='geom' )
    return volclips


def initialise_subsets_in_frame(frame):
    """Function to init all subsets in given frame"""
    poly=fc.get_polygon_from_frame(frame)
    volcs = get_volcanoes_in_polygon(poly)['volc_id'].values
    for v in volcs:
        vids = get_volclip_vids(v)
        for vid in vids:
            initialise_subset_volclip(vid, frame = frame)
    



def initialise_subset_volclip(vid, frame = None, resol_m = 30, full_overlap=True):
    """This will initialise the volclip subset for given frame and volclip.
    If frame is 'None', it will do this for all relevant (fully overlapping) frames.
    The 'full_overlap' param means, only frames that have full overlap with the volclip are to be used. it has no effect if the frame id is given.
    """
    vclipdb = get_volclips_gpd(vid)
    #volcid = vclipdb.volc_id.values[0]
    lon1, lat1, lon2, lat2 = list(vclipdb.geom.bounds.values[0])
    if not frame:
        print('getting related frames')
        frames = fc.subset_get_frames(lon1, lon2, lat1, lat2, full_overlap=full_overlap, only_initialised=True)
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

def get_volcanoes_within(lon, lat, radius_km = 10, volcs = None):
    ''' this will get all volcanoes in radius from given lon/lat'''
    # radius_m = int(round(radius_km*1000))
    # this unfortunately does not work:
    # sql = 'SELECT * FROM volcanoes WHERE ST_DWithin(geometry, POINT({0},{1}), {2});'.format(str(lon), str(lat), str(radius_m))
    radius_deg = radius_km/111.111
    poly = fc.lonlat_to_poly(lon-radius_deg, lon+radius_deg, lat-radius_deg, lat+radius_deg)
    selvolcs = get_volcanoes_in_polygon(poly, volcs)
    return selvolcs


def oneoff_import_volc_portal_names_to_volcdb(volcs = None):
    ''' one-off import of volc portal names to the database...
    should work with limited selection of volcs in pd table as from get_volc_info()'''
    print('Importing volc portal names to the database')
    volcprocdir = '/gws/pw/j07/comet_lics/LiCSAR_volc/volc-proc'
    volctxts = glob.glob(volcprocdir+'/list_database/*volcano*txt')
    if type(volcs)==type(None):
        volcs = get_volc_info()
    ch = []
    for txt in volctxts:
        txtpd = pd.read_csv(txt, header=None, delim_whitespace=True)
        vportal_area = os.path.basename(txt).split('volcano')[0][:-1]
        print(vportal_area)
        #
        for i,vrow in txtpd.iterrows():
            vportal_name = vrow[0]
            vlat = vrow[1]
            vlon = vrow[2]
            selvolc = get_volcanoes_within(vlon, vlat, radius_km=1, volcs=volcs)
            if not len(selvolc) == 0:
                # trying increase range
                selvolc = get_volcanoes_within(vlon, vlat, radius_km=8, volcs=volcs)
            if not len(selvolc)==1:
                ch.append((vportal_name, vlat, vlon))
                print('please check: '+vportal_name+' - found records: '+str(len(selvolc)))
            else:
                volcid=selvolc['volc_id'].values[0]
                sql_q = "UPDATE volcanoes SET vportal_area='{0}',vportal_name='{1}'  where volc_id = {2};".format(vportal_area, vportal_name, str(volcid))
                res = do_query(sql_q, True)
                if not res == 1:
                    print('Error, please check on '+str(volcid))
                    ch.append((vportal_name, vlat, vlon))
    return ch

'''
try then:
import subprocess as sub
volcs=get_volc_info()
volcprocdir = '/gws/pw/j07/comet_lics/LiCSAR_volc/volc-proc'
volctxts = glob.glob(volcprocdir+'/list_database/*volcano*txt')
chs2 = []
for h in chs:
vportal_name = h[0]
vlat, vlon = h[1], h[2]
a = get_volcanoes_within(vlon, vlat, radius_km=9, volcs=volcs)
    if len(a)==1:
        print(vportal_name+' == '+a['name'].values[0])
        vportal_area=sub.getoutput('grep -l '+vportal_name+' '+volcprocdir+'/list_database/*volcano*txt')
        vportal_area=os.path.basename(vportal_area).split('volcano')[0][:-1]
        volcid=a['volc_id'].values[0]
        sql_q = "UPDATE volcanoes SET vportal_area='{0}',vportal_name='{1}'  where volc_id = {2};".format(vportal_area, vportal_name, str(volcid))
        res = do_query(sql_q, True)
    else:
        print(h[0]+' found in '+str(len(a))+' records')
        chs2.append((vportal_name, vlat, vlon))
'''