#!/usr/bin/env python
import os
import LiCSquery as lq
import datetime as dt
import fiona
import geopandas as gpd
from shapely.wkt import loads
from shapely.geometry import Polygon
import matplotlib.pyplot as plt
import LiCSAR_lib.LiCSAR_misc as misc
gpd.io.file.fiona.drvsupport.supported_drivers['KML'] = 'rw'

def vis_aoi(aoi):
    # to visualize a polygon element ('aoi')
    crs = {'init': 'epsg:4326'}
    if type(aoi)==list:
        aoi_gpd = gpd.GeoDataFrame(crs=crs, geometry=aoi)
        #for a in aoi:
    else:
        aoi_gpd = gpd.GeoDataFrame(index=[0], crs=crs, geometry=[aoi])
    # load world borders for background
    world = gpd.read_file(gpd.datasets.get_path('naturalearth_lowres'))
    base = world.plot(color='lightgrey', edgecolor='white')
    aoi_gpd.plot(ax=base, color='None', edgecolor='black')
    bounds = aoi_gpd.geometry.bounds
    plt.xlim([bounds.minx.min()-2, bounds.maxx.max()+2])
    plt.ylim([bounds.miny.min()-2, bounds.maxy.max()+2])
    plt.grid(color='grey', linestyle='-', linewidth=0.2)

def vis_bidtanxs(bidtanxs):
    tovis = []
    for bid in bidtanxs:
        tovis.append(lq.get_polygon_from_bidtanx(bid))
    vis_aoi(tovis)

def vis_frame(frame):
    ai = lq.get_bursts_in_frame(frame)
    bidtanxs = lq.sqlout2list(ai)
    vis_bidtanxs(bidtanxs)

def extract_bursts_by_track(bidtanxs, track):
    newbids = []
    for bidtanx in bidtanxs:
        if bidtanx.split('_')[0] == str(track):
            newbids.append(bidtanx)
    return newbids

def bursts2geopandas(bidtanxs, merge = False):
    # in order to export to KML:
    # frame_gpd.to_file('~/kmls/'+frame+'.kml', driver='KML')
    # or to SHP:
    # frame_gpd.to_file('~/shps/'+frame+'.shp', driver='ESRI Shapefile')
    geometry = []
    crs = {'init': 'epsg:4326'}
    if merge == False:
        if type(bidtanxs)==list:
            for bid in bidtanxs:
                geometry.append(lq.get_polygon_from_bidtanx(bid))
            df_name = {'burstID': bidtanxs}
            aoi_gpd = gpd.GeoDataFrame(df_name, crs=crs, geometry=geometry)
        else:
            geometry.append(lq.get_polygon_from_bidtanx(bidtanxs))
            aoi_gpd = gpd.GeoDataFrame(index=[0], crs=crs, geometry=[geometry])
    else:
        orbdir = lq.get_orbit_from_bidtanx(bidtanxs[0])
        polygon = generate_frame_polygon(bidtanxs, orbdir)
        framename = generate_frame_name(bidtanxs)
        #aoi_gpd = gpd.GeoDataFrame(index=[0], crs=crs, geometry=[polygon])
        aoi_gpd = gpd.GeoDataFrame({'frameID': [framename]}, crs=crs, geometry=[polygon])
    return aoi_gpd

def frame2geopandas(frame):
    bidtanxs = lq.get_bursts_in_frame(frame)
    bidtanxs = lq.sqlout2list(bidtanxs)
    try:
        newname = generate_frame_name(bidtanxs)
    except:
        print('some problem generating frame name from the bursts of frame: '+frame)
        return None
    if frame[-6:] != newname[-6:]:
        #print('WARNING! This frame changed its definition')
        #print('{0} ---> {1}'.format(frame,newname))
        print('framecare_rename.sh {0} {1}'.format(frame,newname))
    gpan = bursts2geopandas(bidtanxs, merge = True)
    #outgpd = outgpd.append(gpan, ignore_index=True)
    return gpan

def export_geopandas_to_kml(gpan, outfile):
    gpan.to_file(outfile, driver='KML')

def bursts_group_to_iws(bidtanxs):
    iw1s = []
    iw2s = []
    iw3s = []
    for bidt in bidtanxs:
        if 'IW1' in bidt:
            iw1s.append(bidt)
        if 'IW2' in bidt:
            iw2s.append(bidt)
        if 'IW3' in bidt:
            iw3s.append(bidt)
    iw1s.sort()
    iw2s.sort()
    iw3s.sort()
    return [iw1s, iw2s, iw3s]

def generate_frame_polygon(bidtanxs, orbdir):
    [iw1s, iw2s, iw3s] = bursts_group_to_iws(bidtanxs)
    if orbdir == 'A':
        #print('not yet tested')
        iwsmin = 0
        iwsmax = -1
        first_point_id = 0
        second_point_id = 1
        last_point = -1
        prelast_point = -2
    else:
        iwsmin = -1
        iwsmax = 0
        first_point_id = 0
        second_point_id = 1
        last_point = -1
        prelast_point = -2
    minbids = []
    maxbids = []
    for iws in [iw1s, iw2s, iw3s]:
        if len(iws)>0:
            minbids.append(iws[iwsmin])
            maxbids.append(iws[iwsmax])
    lons_poly=[]
    lats_poly=[]
    for bid in minbids:
        bidpoly = lq.get_polygon_from_bidtanx(bid)
        xy = bidpoly.exterior.coords.xy
        x = set()
        y = set()
        for pom in xy[0]: x.add(pom)
        for pom in xy[1]: y.add(pom)
        x = list(x)
        y = list(y)
        x.sort()
        y.sort()
            #get two minimal points of lats
        lats_poly.append(y[first_point_id])
        lats_poly.append(y[second_point_id])
        index_min0 = xy[1].index(y[first_point_id])
        index_min1 = xy[1].index(y[second_point_id])
            #get lons that correspond to their lats
        lons_poly.append(xy[0][index_min0])
        lons_poly.append(xy[0][index_min1])
    maxbids.reverse()
    for bid in maxbids:
        bidpoly = lq.get_polygon_from_bidtanx(bid)
        xy = bidpoly.exterior.coords.xy
        x = []
        y = []
        for pom in xy[0]: x.append(pom)
        for pom in xy[1]: y.append(pom)
        x.sort()
        y.sort()
        #get two maximal points of lats
        lats_poly.append(y[last_point])
        lats_poly.append(y[prelast_point])
        index_max0 = xy[1].index(y[last_point])
        index_max1 = xy[1].index(y[prelast_point])
        lons_poly.append(xy[0][index_max0])
        lons_poly.append(xy[0][index_max1])
    return Polygon(zip(lons_poly, lats_poly))

def generate_frame_name(bidtanxs):
    track = bidtanxs[0].split('_')[0]
    orbdir = lq.get_orbit_from_bidtanx(bidtanxs[0])
    polyhon = generate_frame_polygon(bidtanxs, orbdir)
    lat_center = polyhon.centroid.xy[1][0]
    colat = misc.get_colat10(lat_center)
    #print(colat)
    polyid_track = '00'+str(track)+orbdir
    polyid_track = polyid_track[-4:]
    [iw1s, iw2s, iw3s] = bursts_group_to_iws(bidtanxs)
    iw1_str = str(len(iw1s)); iw2_str = str(len(iw2s)); iw3_str = str(len(iw3s))
    if len(iw1s) < 10: iw1_str = '0'+str(len(iw1s))
    if len(iw2s) < 10: iw2_str = '0'+str(len(iw2s))
    if len(iw3s) < 10: iw3_str = '0'+str(len(iw3s))
    polyid_colat10 = str(colat)
    if colat < 10000: polyid_colat10 = '0'+str(colat)
    if colat < 1000: polyid_colat10 = '00'+str(colat)
    if colat < 100: polyid_colat10 = '000'+str(colat)
    if colat < 10: polyid_colat10 = '0000'+str(colat)
    polyid_name = polyid_track+'_'+polyid_colat10+'_'+iw1_str+iw2_str+iw3_str
    return polyid_name

def generate_new_frame(bidtanxs,testonly = True):
    #and now i can generate the new frame:
    track = bidtanxs[0].split('_')[0]
    orbdir = lq.get_orbit_from_bidtanx(bidtanxs[0])
    polyhon = generate_frame_polygon(bidtanxs, orbdir)
    lat_center = polyhon.centroid.xy[1][0]
    colat = misc.get_colat10(lat_center)
    #print(colat)
    polyid_track = '00'+str(track)+orbdir
    polyid_track = polyid_track[-4:]
    [iw1s, iw2s, iw3s] = bursts_group_to_iws(bidtanxs)
    iw1_str = str(len(iw1s)); iw2_str = str(len(iw2s)); iw3_str = str(len(iw3s))
    if len(iw1s) < 10: iw1_str = '0'+str(len(iw1s))
    if len(iw2s) < 10: iw2_str = '0'+str(len(iw2s))
    if len(iw3s) < 10: iw3_str = '0'+str(len(iw3s))
    polyid_colat10 = str(colat)
    if colat < 10000: polyid_colat10 = '0'+str(colat)
    if colat < 1000: polyid_colat10 = '00'+str(colat)
    if colat < 100: polyid_colat10 = '000'+str(colat)
    if colat < 10: polyid_colat10 = '0000'+str(colat)
    polyid_name = polyid_track+'_'+polyid_colat10+'_'+iw1_str+iw2_str+iw3_str
    #print(polyid_name)
    sql = "select count(*) from polygs where polyid_name = '{0}';".format(polyid_name)
    polyid_exists = lq.do_query(sql)[0][0]
    if polyid_exists:
        print('the polyid_name '+ polyid_name +' exists, skipping')
        return False
    lats = []
    lons = []
    for i in range(12):
        lats.append('NULL')
        lons.append('NULL')
    #print(len(polyhon.exterior.coords.xy[1]))
    # removing last coordinate as this should be same as the first one
    for i in range(len(polyhon.exterior.coords.xy[1])-1):
        #lats.append(lat)
        lats[i] = polyhon.exterior.coords.xy[1][i]
        lons[i] = polyhon.exterior.coords.xy[0][i]
    #for lon in polyhon.exterior.coords.xy[0]:
    #    lons.append(lon)
    sql = 'select polyid from polygs order by polyid desc limit 1;'
    lastpolyid = lq.do_query(sql)
    polyid = int(lastpolyid[0][0])+1
    inserted = str(dt.datetime.now())
    name_old = 'frame_'+polyid_track[-1]+'_t'+str(int(polyid_track[0:3]))+'_1bidxxxxxx'
    sql = "INSERT INTO polygs VALUES ({0}, '{1}', '{2}', {3}, {4}, {5}, {6}, {7}, {8}, {9}, {10}, {11}, "\
    "{12}, {13}, {14}, {15}, {16}, {17}, {18}, {19}, {20}, {21}, {22}, {23}, {24}, {25}, {26}, {27}, "\
    "{28}, {29}, {30}, '{31}', {32}, '{33}');".format(\
                       polyid, polyid_name, polyid_track, polyid_colat10, len(iw1s), len(iw2s), len(iw3s),\
                       lats[0], lons[0], lats[1], lons[1], lats[2], lons[2], lats[3], lons[3],\
                       lats[4], lons[4], lats[5], lons[5], lats[6], lons[6], lats[7], lons[7],\
                       lats[8], lons[8], lats[9], lons[9], lats[10], lons[10], lats[11], lons[11],\
                       name_old, 0, inserted)
    #print(sql)
    #return 0
    if testonly:
        print('TEST ONLY')
        print('this command would be performed: ')
        print(sql)
    else:
        res = lq.do_query(sql, 1)
    bid_tanx = iw1s+iw2s+iw3s
    for bidtanx in bid_tanx:
        sql = 'select bid from bursts where bid_tanx="{}";'.format(bidtanx)
        res = lq.do_query(sql)
        #print(sql)
        bid = res[0][0]
        sql = 'INSERT INTO polygs2bursts VALUES ({0}, {1});'.format(polyid,bid)
        #print(sql)
        if testonly:
            print(sql)
        else:
            res = lq.do_query(sql, 1)
    if not testonly:
        print('generated new frame '+polyid_name)
        print('you may do following now: ')
        print('licsar_initiate_new_frame.sh '+polyid_name)
        #delete_frame_commands(frame)
    return polyid_name

def load_bursts_from_txt(intxt):
    f = open(intxt,'r')
    contents = f.readlines()
    bursts = []
    for burst in contents:
        bursts.append(burst.split('\n')[0])
    f.close()
    return bursts

def load_bursts_from_kml(inputkml):
    #inputkml = '/gws/nopw/j04/nceo_geohazards_vol2/LiCS/temp/insar_temp/frames_redef/kmls/test_146a.kml'
    newbursts = gpd.read_file(inputkml, driver='KML')
    newbursts = newbursts[newbursts.columns[0]].tolist()
    return newbursts

def export_bidtanxs_to_kml(bidtanxs, outpath = '/gws/nopw/j04/nceo_geohazards_vol2/LiCS/temp/insar_temp/frames_redef/kmls', projname = 'track', merge = False):
    #kmlout name will be auto_completed
    bidtanxs.sort()
    tracks = set()
    for bid in bidtanxs:
        tracks.add(bid.split('_')[0])
    #print(tracks)
    for track in tracks:
        track_bursts = []
        for bidtanx in bidtanxs:
            if bidtanx.split('_')[0] == track:
                track_bursts.append(bidtanx)
        #print(track_bursts)
        if len(tracks) > 1:
            kmlout = os.path.join(outpath,'{0}_{1}.kml'.format(projname, track))
        else:
            kmlout = os.path.join(outpath,'{0}.kml'.format(projname))
        #print(kmlout)
        if os.path.exists(kmlout): os.remove(kmlout)
        frame_gpd = bursts2geopandas(track_bursts, merge)
        print('exporting to '+kmlout)
        frame_gpd.to_file(kmlout, driver='KML')
    if merge == False:
        print('done. please edit the kmls - delete not wanted bursts, save and return')

def export_frame_to_kml(frame, outpath = '/gws/nopw/j04/nceo_geohazards_vol2/LiCS/temp/insar_temp/frames_redef/kmls'):
    if not os.path.exists(outpath):
        os.mkdir(outpath)
    ai = lq.get_bursts_in_frame(frame)
    bidtanxs = lq.sqlout2list(ai)
    export_bidtanxs_to_kml(bidtanxs, outpath, projname = frame)

def delete_frame_commands(frame):
    print('setFrameInactive.py {0}'.format(frame))
    print('rm -rf $LiCSAR_procdir/{0}/{1} $LiCSAR_public/{0}/{1}'.format(track,frame))
    print('lics_mysql.sh')
    sql = "select polyid from polygs where polyid_name = '{0}';".format(frame)
    polyid = lq.do_query(sql)[0][0]
    sql = "delete from polygs2master where polyid={};".format(polyid)
    print(sql)
    sql = "delete from polygs2bursts where polyid={};".format(polyid)
    print(sql)
    frame_workaround = frame.replace('A','X')
    frame_workaround = frame_workaround.replace('D','X')
    sql = "update polygs set polyid_name='{0}' where polyid_name='{1}';".format(frame_workaround, frame)
    print(sql)
    sql = "delete from polygs where polygs.polyid_name='{}';".format(frame_workaround)
    print(sql)
    print('')
    track=str(int(frame[0:3]))

