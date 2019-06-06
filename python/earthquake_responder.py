import LiCSquery as lq
import datetime
from shapely.geometry import Polygon
from shapely.ops import unary_union

def get_bursts_within_polygon(lon1,lat1,lon2,lat2):
    sql_q = "select bid, bid_tanx, bid_tanxtk "\
            'corner1_lon, corner1_lat, '\
            'corner2_lon, corner2_lat, '\
            'corner3_lon, corner3_lat, '\
            'corner4_lon, corner4_lat '\
            "from bursts " \
            "where "\
            "greatest( "\
            "corner1_lon,corner2_lon,corner3_lon,corner4_lon"\
            ") "\
            ">= {0} and "\
            "least( "\
            "corner1_lon,corner2_lon,corner3_lon,corner4_lon"\
            ") "\
            "<= {1} ".format(lon1,lon2)
    sql_q += "and "\
            "greatest( "\
            "corner1_lat,corner2_lat,corner3_lat,corner4_lat"\
            ") "\
            ">= {0} and "\
            "least( "\
            "corner1_lat,corner2_lat,corner3_lat,corner4_lat"\
            ") "\
            "<= {1};".format(lat1,lat2)
    #print(sql_q)
    return lq.do_query(sql_q)

def get_range_from_magnitude(mag,depth):
    #from alistair codes:
    #but it is weird... too small value
    #on the other hand, John's table showed 1150 km - and this is way too much ..
    #...or not - USGS showed that IV shaking was felt some 750 km from the epicenter!
    ang = 10 + mag * 10
    radius = depth * np.tan(ang/2)
    return radius

def get_polygon_from_burst_id(bid):
    sql = 'select corner1_lat, corner1_lon, corner2_lat, corner2_lon, corner3_lat, corner3_lon,'\
    ' corner4_lat, corner4_lon from bursts where bid = {0};'.format(bid)
    coords = lq.do_query(sql)[0]
    #the coordinates are 'random' so I do a loop to get valid polygon
    lat_set = [(0, 2, 4, 6),(0, 2, 6, 4),(0, 4, 6, 2),(0, 4, 2, 6),(0, 6, 4, 2),(0, 6, 2, 4)]
    lon_set = [(1, 3, 5, 7),(1, 3, 7, 5),(1, 5, 7, 3),(1, 5, 3, 7),(1, 7, 5, 3),(1, 7, 3, 5)]
    for i in range(len(lat_set)):
        lat_point_list = []
        lon_point_list = []
        for x in lat_set[i]:
            lat_point_list.append(coords[x])
        for y in lon_set[i]:
            lon_point_list.append(coords[y])
        polygon_geom = Polygon(zip(lon_point_list, lat_point_list))
        if polygon_geom.is_valid:
            break
    return polygon_geom

def create_new_frame(burst_ids,iw1,iw2,iw3,polyid_track,eq_date):
    #if polyid_name exists in db, just return its name. otherwise:
    #insert into polygs: polyid_name=011D_00000_050505
    datestring = eq_date.strftime("%y%m%d")[1:]
    iw1_str = str(iw1); iw2_str = str(iw2); iw3_str = str(iw3)
    if iw1 < 10: iw1_str = '0'+str(iw1)
    if iw2 < 10: iw2_str = '0'+str(iw2)
    if iw3 < 10: iw3_str = '0'+str(iw3)
    if (iw1 > 99) or (iw2 > 99) or (iw3 > 99):
        print('The track '+polyid_track+' exceeds limit of less than 100 bursts per swath')
        return
    polyid_name = polyid_track+'_'+datestring+'_'+iw1_str+iw2_str+iw3_str
    #check if such polyid_name exists:
    sql = "select count(*) from polygs where polyid_name = '{0}';".format(polyid_name)
    polyid_exists = lq.do_query(sql)[0][0]
    if polyid_exists:
        print('the polyid_name '+ polyid_name +' exists, skipping')
        return
    polyid_colat10 = datestring
    name_old = 'eqake_'+polyid_track[-1]+'_t'+str(int(polyid_track[0:3]))+'_'+eq_date.strftime("%y%m%d")
    #hardest part: to identify corner lat, lon coordinates
    print('Getting frame polygon from bursts')
    polygons = []
    for bid in burst_ids:
        polygons.append(get_polygon_from_burst_id(bid))
    merged_polygon = unary_union(polygons)
    #this workaround helps if we have some gaps between bursts
    print(polyid_name)
    if merged_polygon.geom_type == 'MultiPolygon':
        merged_polygon = merged_polygon.convex_hull
    #print('number of original points: '+str(len(merged_polygon.exterior.coords)-2))
    for simp in range(1,10):
        frame_polygon = merged_polygon.simplify(0.1*simp, preserve_topology=True)
        #print(str(len(frame_polygon.exterior.coords)-2))
        if len(frame_polygon.exterior.coords)-2 <= 12:
            break
    #if still complex, do hull and simplify again
    if len(frame_polygon.exterior.coords)-2 > 12:
        frame_polygon = frame_polygon.convex_hull
    if len(frame_polygon.exterior.coords)-2 > 12:
        for simp in range(1,100):
            frame_polygon = frame_polygon.simplify(0.1*simp, preserve_topology=False)
            if len(frame_polygon.exterior.coords)-2 <= 12:
                break
    #boundary = gpd.GeoSeries(merged_polygon)
    #boundary.plot(color = 'blue')
    #plt.show()
    #print('number of simplified points: '+str(len(frame_polygon.exterior.coords)-2))
    #boundary = gpd.GeoSeries(frame_polygon)
    #boundary.plot(color = 'red')
    #plt.show()
    #return
    max_i = len(frame_polygon.exterior.coords)-1
    #get minimum number of 4 vertices
    if max_i < 4: max_i = 4
    lons = []; lats = []
    for i in range(12):
        lons.append('NULL')
        lats.append('NULL')
    for i in range(1, max_i):
        lons[i-1] = frame_polygon.exterior.coords[i][0]
        lats[i-1] = frame_polygon.exterior.coords[i][1]
    #now to insert the new temporary frame to the database
    sql = 'select polyid from polygs order by polyid desc limit 1;'
    lastpolyid = lq.do_query(sql)
    polyid = int(lastpolyid[0][0])+1
    inserted = str(datetime.datetime.now())
    sql = "INSERT INTO polygs VALUES ({0}, '{1}', '{2}', {3}, {4}, {5}, {6}, {7}, {8}, {9}, {10}, {11}, "\
    "{12}, {13}, {14}, {15}, {16}, {17}, {18}, {19}, {20}, {21}, {22}, {23}, {24}, {25}, {26}, {27}, "\
    "{28}, {29}, {30}, '{31}', {32}, '{33}');".format(\
                       polyid, polyid_name, polyid_track, polyid_colat10, iw1, iw2, iw3,\
                       lats[0], lons[0], lats[1], lons[1], lats[2], lons[2], lats[3], lons[3],\
                       lats[4], lons[4], lats[5], lons[5], lats[6], lons[6], lats[7], lons[7],\
                       lats[8], lons[8], lats[9], lons[9], lats[10], lons[10], lats[11], lons[11],\
                       name_old, 0, inserted)
    #print(sql)
    res = lq.do_query(sql, 1)
    #and fill the polygs2bursts for relation polyid<->burst_ids
    for bid in burst_ids:
        sql = 'INSERT INTO polygs2bursts VALUES ({0}, {1});'.format(polyid,bid)
        #print(sql)
        res = lq.do_query(sql, 1)
    return polyid_name





#and this is the main code:
#... so this code actually defines polygons (frames) containing the eathquake
# i use now only licsar_live ...
# probably should change it - and clean the polygons that have polyid > 9000

#center coordinates:
#let's try the Peru earthquake:
lat = -5.807
lon = -75.264
MAG = 8
depth = 122 #km
rad_km = get_range_from_magnitude(MAG,depth)
#range in km:
rad_km = 750
#date of earthquake:
#to make things go properly, we should introduce also time here!
eq_date = datetime.datetime(2019,5,26)
radius = (360.0 / 40007.86) * rad_km
minlon = lon - radius
maxlat = lat + radius
maxlon = lon + radius
minlat = lat - radius
selected_bursts = get_bursts_within_polygon(minlon,minlat,maxlon,maxlat)
#print(selected_bursts)
print(str(len(selected_bursts))+' bursts selected')

burst_ids=[]
relorbs=[]
for i in range(len(selected_bursts)):
    burst_ids.append(selected_bursts[i][1])
    relorbs.append(selected_bursts[i][2])
relorbs = list(dict.fromkeys(relorbs))
print(str(len(relorbs))+' orbital tracks')

polyids = []
#one of polyids as i tested was: 127D_90526_285070
for relorb in relorbs:
    bids=[]
    iw1=0
    iw2=0
    iw3=0
    pom=0
    #datestring = eq_date.strftime("%y%m%d")[1:]
    for i in range(len(selected_bursts)):
        if selected_bursts[i][2] == relorb:
            if pom==0:
                #i should get orb_dir here
                bid=str(selected_bursts[i][0])
                sql='select f.orb_dir from files2bursts fb inner join files f on fb.fid=f.fid where fb.bid={0} limit 1;'.format(bid)
                orb_dir=lq.do_query(sql)[0][0]
                relorb_str = str(relorb)
                if relorb < 100: relorb_str = '0'+relorb_str
                if relorb < 10: relorb_str = '0'+relorb_str
                polyid_track = relorb_str+orb_dir
                #print(orb_dir)
                pom=1
            bids.append(selected_bursts[i][0])
            iw=int(selected_bursts[i][1].split('_')[1][2])
            if iw==1: iw1+=1
            if iw==2: iw2+=1
            if iw==3: iw3+=1
    #print(bids)
    #print(iw1)
    polyids.append(create_new_frame(bids,iw1,iw2,iw3,polyid_track,eq_date))
print(polyids)
