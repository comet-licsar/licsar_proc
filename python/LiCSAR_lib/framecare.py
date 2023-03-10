#!/usr/bin/env python

# Mar 2020+ - Milan Lazecky

import os, glob
import subprocess as subp
import LiCSquery as lq
from LiCSAR_misc import *
import datetime as dt
import fiona
import pandas as pd
import geopandas as gpd
from shapely.wkt import loads
from shapely.geometry import Polygon
import matplotlib.pyplot as plt
import LiCSAR_lib.LiCSAR_misc as misc
import s1data as s1
gpd.io.file.fiona.drvsupport.supported_drivers['KML'] = 'rw'
import rioxarray

pubdir = os.environ['LiCSAR_public']
procdir = os.environ['LiCSAR_procdir']


'''
# notes:
# this is how i imported burst db - first i converted them from sqlite3 to geojson
# then i did, in /gws/nopw/j04/nceo_geohazards_vol1/projects/LiCS/proc/current/burst_database/IW/sqlite:

import geopandas as gpd
import shapely
import time
import LiCSquery as lq

aa=gpd.read_file('burst_map.geojson')
aa=aa[aa.burst_id>56099]

aa['geometry']=aa['geometry'].convex_hull

def _to_2d(x, y, z):
    return tuple(filter(None, [x, y]))

aa['geometry'] = aa['geometry'].apply(lambda x: shapely.ops.transform(_to_2d, x))

# now it is ready to import to database:
for i,j in aa.iterrows():
    print(i)
    res = store_burst_geom(j[0], int(j[1][-1]), j[2], j[3], j[4][0], j[5].wkt)
    #time.sleep(0.25)

lon1=72.510
lon2=72.845
lat1=38.130
lat2=38.365
resol_m=30
frame='100A_05236_141313'
sid='SAREZ'
fc.subset_initialise_corners(frame, lon1, lon2, lat1, lat2, sid)


'''

def subset_initialise_corners(frame, lon1, lon2, lat1, lat2, sid, is_volc = False, resol_m=30):
    """This will initialise a subset given by corner lon/lat-s.
    The results will be stored in $LiCSAR_procdir/subsets
    
    Args:
        frame (str): frame ID,
        lon1, lon2 (float, float): corner longitudes,
        lat1, lat2 (float, float): corner latitudes,
        sid (str):  string ID (for volcano, use its volcano ID number)
        is_volc (bool): if true, it will set the output folder $LiCSAR_procdir/subsets/volc
        resol_m (float): output resolution in metres to have geocoding table ready in (note, RSLCs are anyway in full res)
    """
    if is_volc:
        sidpath = 'volc/'+sid
    else:
        sidpath = sid
    #
    resol=resol_m/111111 #0.00027
    resol=round(resol,6)
    # sort the coordinates
    lon1,lon2=sorted([lon1,lon2])
    lat1,lat2=sorted([lat1,lat2])
    #
    track=str(int(frame[0:3]))
    framedir = os.path.join(os.environ['LiCSAR_procdir'],track,frame)
    subsetdir = os.path.join(os.environ['LiCSAR_procdir'],'subsets',sidpath,frame[:4])
    if os.path.exists(subsetdir):
        print('the subset directory exists. continuing anyway..')
    if not os.path.exists(os.path.join(framedir, 'subsets')):
        os.mkdir(os.path.join(framedir, 'subsets'))
    #
    # get median height
    print('getting median height')
    hgt=os.path.join(os.environ['LiCSAR_public'], str(int(frame[:3])), frame, 'metadata', frame+'.geo.hgt.tif')
    a=rioxarray.open_rasterio(hgt)
    a=a.sortby(['x','y'])
    medhgt=round(float(a.sel(x=slice(lon1,lon2), y=slice(lat1, lat2)).median()))
    #medhgt=round(float(a.sel(x=(lon1,lon2), y=(lat1, lat2), method='nearest').median()))
    print('... as {} m'.format(str(medhgt)))
    #
    # running the clipping in init-only mode
    clipcmd = "cd "+framedir+"; "
    clipcmd = clipcmd + "clip_slc.sh "+subsetdir+" "+str(lon1)+" "+str(lon2)+" "
    clipcmd =     clipcmd   +str(lat1)+" "+str(lat2)+" "
    clipcmd =     clipcmd   +str(medhgt)+" "+str(resol)+" 0 1"
    #
    print('initializing the subset')
    os.chdir(framedir)
    print(clipcmd)
    os.system(clipcmd)
    if os.path.exists(subsetdir):
        subsetlink = os.path.join(framedir, 'subsets', sid)
        if not os.path.exists(subsetlink):
            os.symlink(subsetdir, subsetlink)
    else:
        print('some error occurred and the output dir was not created')
    return


def subset_initialise_centre_coords(frame, clon, clat, sid, is_volc = False, radius = 25/2, resol_m=30):
    """This will initialise a subset given by centre lon/lat and radius in km.
    The results will be stored in \$LiCSAR_procdir/subsets
    
    Args:
        frame (str): frame ID,
        clon (float): centre longitude,
        clat (float): centre latitude,
        sid (str):  string ID (for volcano, use its volcano ID number)
        is_volc (bool): if true, it will set the output folder \$LiCSAR_procdir/subsets/volc
        radius (float): radius (half of the diameter) of the subset scene, in km
        resol_m (float): output resolution in metres to have geocoding table ready in (note, RSLCs are anyway in full res)
    """
    
    # BASICS
    radius_deg=radius_km/111
    lon1=lon-radius_deg
    lon2=lon+radius_deg
    lat1=lat-radius_deg
    lat2=lat+radius_deg
    subset_initialise_corners(frame, lon1, lon2, lat1, lat2, sid, is_volc = False, resol_m=30)
    return


def check_and_fix_burst(mburst, framebursts):
    # to get mbursts of a zip file, e.g.:
    # frame = '...'
    # filename = 'S1A_IW_SLC__1SDV_20210429T114802_20210429T114829_037665_047199_F24F.zip'
    # mbursts = fc.lq.sqlout2list(fc.lq.get_bursts_in_file(filename))
    # framebursts = fc.lq.sqlout2list(fc.lq.get_bidtanxs_in_frame(frame))
    # for mburst in mbursts: fc.check_and_fix_burst(mburst, framebursts)
    changed = False
    if mburst in framebursts:
        return changed
    tr=int(mburst.split('_')[0])
    iw=mburst.split('_')[1]
    tanx=int(mburst.split('_')[2])
    #
    for fburst in framebursts:
        iwf=fburst.split('_')[1]
        if iwf == iw:
            tanxf=int(fburst.split('_')[2])
            # checking in a 'relaxed' tolerance (0.8 s)
            if abs(tanx - tanxf) < 8:
                # just to make sure they are both of the same pass..
                if lq.get_orbdir_from_bidtanx(fburst) == lq.get_orbdir_from_bidtanx(mburst):
                    # check if their geometries overlap
                    fb_gpd = bursts2geopandas([fburst])
                    mb_gpd = bursts2geopandas([mburst])
                    if fb_gpd.overlaps(mb_gpd).values[0]:
                        print('we (very) probably found a cross-defined burst. fixing/merging to one')
                        print(mburst+' -> '+fburst)
                        lq.rename_burst(mburst, fburst)
                        changed = True
    return changed



def check_and_fix_all_bursts_in_frame(frame):
    t1 = '2014-10-01'
    t2 = dt.datetime.now().date()
    framefiles = lq.get_frame_files_period(frame,t1,t2)
    framebursts = lq.sqlout2list(lq.get_bidtanxs_in_frame(frame))
    fdates = []
    noch = 0
    for framefile in framefiles:
        filename=framefile[2]+'.zip'
        print('checking '+filename)
        mbursts = lq.sqlout2list(lq.get_bursts_in_file(filename))
        # visual check
        # from matplotlib import pyplot as plt
        # bursts_gpd = bursts2geopandas(mbursts)
        # frame_gpd = bursts2geopandas(framebursts)
        # bursts_gpd.plot()
        # frame_gpd.plot()
        # plt.show()

        for mburst in mbursts:
            changed = check_and_fix_burst(mburst, framebursts)
            if changed:
                noch = noch + 1
                fdate = filename[17:25]
                fdates.append(fdate)
    fdates = list(set(fdates))
    print('additionally checking only burst ids of similar tracks')
    print('(same orbit pass direction, relorb+-1)')
    #trackid = frame[:4]
    #for fburst in framebursts:
    track = int(frame[:3]) #int(fburst.split('_')[0])
    for cant in [str(track-1), str(track), str(track+1)]:
        if cant == '176':
            cant = '001'
        if cant == '0':
            cant = '175'
        if len(cant) == 1:
            cant = '00'+cant
        if len(cant) == 2:
            cant = '0'+cant
        trackid = cant+frame[3]
        #canfburst = str(cant)+'_'+fburst.split('_')[1]+'_'+fburst.split('_')[2]
        canfbursts = lq.get_bidtanxs_in_track(trackid)
        try:
            canfbursts = lq.sqlout2list(canfbursts)
            for canfburst in canfbursts:
                changed = check_and_fix_burst(canfburst, framebursts)
                if changed:
                    noch = noch + 1
        except:
            print('no bursts in the trackid '+trackid)
    print(str(noch)+' burst definitions changed to fit the frame burst IDs')
    if noch > 0:
        print('you may want to check following epochs:')
        for fdate in fdates:
            print('frame {0}: {1}'.format(frame, fdate))
            #print('remove_from_lics.sh {0} {1}'.format(frame, fdate))

'''
check these frames:
['149D_05278_131313', '150D_05107_131313', '150D_05306_131313', '151D_05241_131313']
['149D_05425_060707', '150D_05306_131313', '150D_05505_131313', '151D_05440_131313']
['149D_05278_131313', '150D_05107_131313', '150D_05306_131313', '151D_05241_131313']

'''

def check_and_fix_all_files_in_frame(frame):
    t1 = '2014-10-01'
    t2 = dt.datetime.now().date()
    files = lq.get_frame_files_period(frame, t1, t2, only_file_title = True)
    files = lq.sqlout2list(files)
    i = 0
    lenf = len(files)
    for f in files:
        i = i+1
        print('['+str(i)+'/'+str(lenf)+'] checking file '+f)
        check_bursts_in_file(f)


def check_and_fix_burst_supershifts_in_frame(frame, viewerror = True):
    frame_wkt = lq.geom_from_polygs2geom(frame)
    framepoly = loads(frame_wkt)
    framepoly_gpd = gpd.GeoSeries(framepoly)
    frame_bursts = lq.sqlout2list(lq.get_bidtanxs_in_frame(frame))
    fgpd = bursts2geopandas(frame_bursts)
    #b1 = fgpd.iloc[0]
    # 4.5 degrees in WGS-84 are approx 500 km - that should be enough to compare from the frame polygon centroid
    cluster1 = fgpd[fgpd.geometry.centroid.distance(framepoly.centroid) <= 4.5]
    cluster2 = fgpd[fgpd.geometry.centroid.distance(framepoly.centroid) > 4.5]
    #polyid = lq.sqlout2list(lq.get_frame_polyid(frame))[0]
    if not cluster2.empty:
        print('here we are - two burst clusters!')
        # checking for the overlap anyways
        if not framepoly.overlaps(cluster1.unary_union):
            badbursts_gpd = cluster1
            goodbursts_gpd = cluster2
        elif not framepoly.overlaps(cluster2.unary_union):
            badbursts_gpd = cluster2
            goodbursts_gpd = cluster1
        if viewerror:
            print('this is how the frame should look like:')
            framepoly_gpd.plot()
            plt.show()
            print('and this is how it looks with the current burst definitions')
            vis_frame(frame)
        #print('trying to solve it - first find one file that has the burst as bad one')
        check_and_fix_all_files_in_frame(frame)
        '''
        for bid in badbursts_gpd.burstID.values:
            filewithbidasbad = ''
            repeat = True
            filescheck = files.copy()
            while repeat:
                filewithbidasbad = ''
                for fileid in filescheck:
                    print(fileid)
                    is_bid_bad_there = check_bursts_in_file(fileid, badburstfind = bid)
                    if is_bid_bad_there == 'yes':
                        #now check if the file has some of the good bids, i.e. if it really is part of the frame
                        filebursts = lq.sqlout2list(lq.get_bursts_in_file(fileid))
                        for fbur in filebursts:
                            if fbur in goodbursts_gpd.burstID.values:
                                filewithbidasbad = fileid
                                break
                        if filewithbidasbad:
                            break
                if not filewithbidasbad:
                    print('no file with this burst as bad one, skipping')
                    repeat = False
                    continue
                #check_bursts_in_file(filewithbidasbad)
                lq.delete_file_from_db(filewithbidasbad, 'name')
                #print('debug')
                #print(filewithbidasbad)
                filepath = s1.get_neodc_path_images(filewithbidasbad, file_or_meta = True)[0]
                #this should regenerate the missing burst
                outchars = ingest_file_to_licsinfo(filepath)
                if outchars:
                    if outchars<200:
                        print('did not help')
                        filescheck.remove(filewithbidasbad)
                        repeat = True
                    else:
                        print('YESSS, the reingested image created new bursts!!!!!')
                        repeat = False
        '''
        # now all the missing bursts probably exist, so let's try getting them in the frame overlap + check colat and exchange bids
        minlon, minlat, maxlon, maxlat = framepoly.bounds
        track = int(frame[:3])
        burstcands = []
        for relorb in [track-1, track, track+1]:
            if relorb == 0: relorb = 175
            if relorb == 176: relorb = 1
            burstcandsT = lq.sqlout2list(lq.get_bursts_in_polygon(minlon, maxlon, minlat, maxlat, relorb))
            for b in burstcandsT:
                burstcands.append(b)
        # now check their number etc.
        #and if all ok, use them instead of the bad ones - replace them
        frame_bursts_to_change = []
        for b in frame_bursts:
            if not b in burstcands:
                frame_bursts_to_change.append(b)
            else:
                burstcands.remove(b)
        if len(frame_bursts_to_change)>len(burstcands):
            print('ERROR - not enough burst candidates - cannot exchange all bursts, cancelling')
            return False
        else:
            #for swath in [1,2,3]:
            frame_bursts_to_change_out = frame_bursts_to_change.copy()
            for fb in frame_bursts_to_change:
                print('checking burst '+fb)
                sw = fb.split('_')[1]
                tanx = fb.split('_')[2]
                for bc in burstcands:
                    if sw == bc.split('_')[1]:
                        if abs(int(bc.split('_')[2])-int(tanx)) < 10:
                            print('exchanging {0} -> {1}'.format(fb,bc))
                            lq.replace_bidtanx_in_frame(frame, fb, bc)
                            frame_bursts_to_change_out.remove(fb)
                            burstcands.remove(bc)
                            break
            if len(frame_bursts_to_change_out) > 0:
                print('ERROR - not all frame bursts were replaced - the problematic bursts are returned:')
                print(frame_bursts_to_change_out)
                print('potential burst candidates were:')
                print(burstcands)
                return [frame_bursts_to_change_out, burstcands]
            else:
                print('the frame was corrected properly!')
                if viewerror:
                    print('see yourself')
                    vis_frame(frame)
            #return burstcands, frame_bursts_to_change
    else:
        print('bursts of this frame are ok')
        return True

# to get all files that are not participating in any frame
#sql = "select f.name from files f where f.fid not in ( select fb.fid from files2bursts fb inner join polygs2bursts pb on fb.bid=pb.bid );"

import time
def process_all_frames():
    badtracks = []
    for relorb in range(1,175): #,175): #   86 need to do: 97-99
        try:
            print('preparing track '+str(relorb+1))
            allframes = lq.sqlout2list(lq.get_frames_in_orbit(relorb+1))
        except:
            print('error in relorb '+str(relorb+1))
            badtracks.append(relorb+1)
            continue
        for frame in allframes:
            print(frame)
            time.sleep(45)
            #just change the function here
            try:
                #rc = check_and_fix_burst_supershifts_in_frame(frame, viewerror = False)
                # to process ALL FILES! (that are related to some any frame)
                rc = check_and_fix_burst_supershifts_in_frame_files(frame, viewerror = False)
            except:
                print('some error during processing frame '+frame)
    return badtracks


'''
#to get files that are NOT in any frames - we have now over 250k of such files!
sql = "select f.name from files f where f.fid not in ( select fb.fid from files2bursts fb inner join polygs2bursts pb on fb.bid=pb.bid );"
nopolyfiles = lq.do_pd_query(sql)
i=0
filez = nopolyfiles.name.unique()
lenn = len(filez)
for fileid in filez:
    i=i+1
    print('['+str(i)+'/'+str(lenn)+']'+fileid)
    lq.delete_file_from_db(fileid, col = 'name')
    filepath = s1.get_neodc_path_images(fileid, file_or_meta = True)[0]
    chars = ingest_file_to_licsinfo(filepath)
    print(chars)
    #time.sleep(2)
    #if not check_bursts_in_file(fileid):
    #    print('error in file '+fileid)
'''

def check_and_fix_burst_supershifts_in_frame_files(frame, viewerror = False, force_reingest = True):
    #first get all files in frame and check them one by one:
    t1 = '2014-10-01'
    t2 = dt.datetime.now().date()
    files = lq.get_frame_files_period(frame, t1, t2, only_file_title = True)
    files = lq.sqlout2list(files)
    for fileid in files:
        if force_reingest:
            print('reingesting '+fileid)
            reingest_file(fileid)
        else:
            if not check_bursts_in_file(fileid):
                print('error in file '+fileid)
                if viewerror:
                    print('see yourself the current situation')
                    bursts = lq.sqlout2list(lq.get_bursts_in_file(fileid))
                    vis_bidtanxs(bursts)


def reingest_all_files_in_frame(frame):
    t1 = '2014-10-01'
    t2 = dt.datetime.now().date()
    files = lq.get_frame_files_period(frame, t1, t2, only_file_title = True)
    files = lq.sqlout2list(files)
    for fileid in files:
        reingest_file(fileid)


def vis_file(fileid):
    """ Visualize bursts of given file.
    """
    fbursts = lq.get_bursts_in_file(fileid)
    fbursts = lq.sqlout2list(fbursts)
    #filegpd = bursts2geopandas(fbursts)
    vis_bidtanxs(fbursts)


def check_bursts_in_file(fileid = 'S1A_IW_SLC__1SDV_20210908T235238_20210908T235305_039597_04AE3C_4CA7', badburstfind = None, autocorrect = True):
    fbursts = lq.get_bursts_in_file(fileid)
    fbursts = lq.sqlout2list(fbursts)
    if not fbursts:
        print('no bursts found for this file. trying to reingest it')
        ingest_file_to_licsinfo(fileid, False)
        return False
    filegpd = bursts2geopandas(fbursts)
    b1 = filegpd.iloc[0]
    # not perfect solution - there can be more than 2 clusters!!!!!!
    # but now, just removing and reingesting the file should help, anyway
    # 4.5 degrees in WGS-84 are approx 500 km - that should be enough..
    cluster1 = filegpd[filegpd.geometry.centroid.distance(b1.geometry.centroid) <= 4.5]
    cluster2 = filegpd[filegpd.geometry.centroid.distance(b1.geometry.centroid) > 4.5]
    if not cluster2.empty:
        print('here we are - two burst clusters!')
        info = s1.get_info_pd(fileid)
        try:
            filepoly = loads(info.footprint.values[0])
        except:
            print('some error loading footprint from scihub')
            print('(making sure things work fine - reingesting this file)')
            filepath = s1.get_neodc_path_images(fileid, file_or_meta = True)[0]
            chars = ingest_file_to_licsinfo(filepath)
            return True
        if not filepoly.overlaps(cluster1.unary_union):
            badbursts_gpd = cluster1
        elif not filepoly.overlaps(cluster2.unary_union):
            badbursts_gpd = cluster2
        else:
            print('weird - both clusters overlap with the original file')
            print('(making sure things work fine - reingesting this file)')
            filepath = s1.get_neodc_path_images(fileid, file_or_meta = True)[0]
            chars = ingest_file_to_licsinfo(filepath)
            return False
        if badburstfind:
            if badburstfind in badbursts_gpd.burstID.values:
                print('this burst is indeed in badbursts')
                return 'yes'
            else:
                return 'no'
        if not autocorrect:
            return badbursts_gpd
        else:
            allremoved = True
            for bid in badbursts_gpd.burstID.values:
                frames = lq.sqlout2list(lq.get_frames_with_burst(bid))
                print('checking '+bid)
                #print(frames)
                if len(frames) == 0:
                    print('no frame is using this burst ID. as it is a bad burst, will remove it now, including files that use it')
                    files2remove = lq.sqlout2list(lq.get_filenames_from_burst(bid))
                    print('removing and reingesting file {}'.format(str(len(files2remove))))
                    for ff in files2remove:
                        lq.delete_file_from_db(ff, col = 'name')
                        filepath = s1.get_neodc_path_images(ff, file_or_meta = True)[0]
                        chars = ingest_file_to_licsinfo(filepath)
                    print('removing the bad burst '+bid)
                    rc = lq.delete_burst_from_db(bid)
                else:
                    #print('this burst is used in following frame(s):')
                    #print(frames)
                    allremoved = False
                    # ye.. or better append and remove dup. but who cares..
                    badframes = frames
            if allremoved:
                print('bad bursts are cleaned! reingesting the file')
            else:
                print('not all bad bursts removed from the database - but reingesting the file to fix it')
                print('SOME of frames still using a bad burst:')
                print(badframes)
            lq.delete_file_from_db(fileid, col = 'name')
            filepath = s1.get_neodc_path_images(fileid, file_or_meta = True)[0]
            chars = ingest_file_to_licsinfo(filepath)
            return True
    else:
        print('bursts of this file are ok')
        return True

'''
filez=filez[4514:]
lenn = len(filez)
i=0
badones = []
for fileid in filez.name.values:
    i=i+1
    print('['+str(i)+'/'+str(lenn)+'] '+fileid)
    if 'IW' in fileid:
        try:
            reingest_file(fileid)
        except:
            print('error with '+fileid)
            badones.append(fileid)
        
        
    time.sleep(5)
    if not check_bursts_in_file(fileid):
        print('error in file '+fileid)

'''

def reingest_file(fileid):
    rc = lq.delete_file_from_db(fileid, col = 'name')
    chars = ingest_file_to_licsinfo(fileid, False)
    return chars


def ingest_file_to_licsinfo(filepath, isfullpath = True):
    """ Will ingest a S1 SLC zip file to the LiCSInfo database.
    If filepath is only filename, it will try find this file in neodc or LiCSAR_SLC"""
    if not isfullpath:
        filepath = s1.get_neodc_path_images(filepath, file_or_meta = True)[0]
    if not os.path.exists(filepath):
        filepath = os.path.join(os.environ['LiCSAR_SLC'], os.path.basename(filepath))
        if not os.path.exists(filepath):
            print('ERROR - this file does not exist')
            return False
        aaaa = subp.check_output(['arch2DB.py','-f',filepath])
    else:
        #cmd = 'arch2DB.py -f {} >/dev/null 2>/dev/null'.format(filepath)
        #cmd = 'arch2DB.py -f {}'.format(filepath)
        #rc = os.system(cmd)
        aaaa = subp.check_output(['arch2DB.py','-f',filepath])
    return len(aaaa)


def get_bidtanxs_from_xy(lon,lat,relorb=None,swath=None, tol=0.05):
    """Gets bursts in given coordinates (and optionally in given track or swath)"""
    bursts = lq.get_bursts_in_xy(lon,lat,relorb,swath,tol)
    bursts = lq.sqlout2list(bursts)
    return bursts


def get_bidtanxs_from_xy_file(intxt, relorb = None):
    """Gets bursts in polygon given by the xy text file."""
    if not os.path.exists(intxt):
        print('ERROR, the file does not exist')
        return False
    lonlat = load_xy(intxt)
    bidtanxs = lq.get_bursts_in_polygon(lonlat[0][0],lonlat[0][-1],lonlat[1][0], lonlat[1][-1], relorb = relorb)
    bidtanxs = lq.sqlout2list(bidtanxs)
    print('check the bursts, e.g. export_bidtanxs_to_kml')
    return bidtanxs


def make_bperp_file(frame, bperp_file):
    """Creates baselines file for given frame, by requesting info from ASF"""
    mid = get_master(frame, asfilenames = True)
    if not mid:
        return False
    mid=mid[0].split('.')[0]
    bpd = s1.get_bperps_asf(mid)
    bpd.to_csv(bperp_file, sep = ' ', index = False, header = False)


def get_master(frame, asfilenames = False, asdate = False, asdatetime = False, metafile = None):
    if not metafile:
        track=str(int(frame[0:3]))
        metafile = os.path.join(pubdir,track,frame,'metadata','metadata.txt')
    if not os.path.exists(metafile):
        print('frame {} is not initialised'.format(frame))
        return False
    master = misc.grep1line('master',metafile)
    if not master:
        print('error parsing information from metadata.txt')
        return False
    masterdate = master.split('=')[1]
    if asfilenames:
        slcpath = os.path.join(procdir, track, frame, 'SLC', str(masterdate))
        try:
            zipfiles = []
            for zipfile in glob.glob(slcpath+'/S1*zip'):
                zipfiles.append(os.path.basename(zipfile))
            return zipfiles
        except:
            print('error finding zip files in the frame SLC directory')
            return False
    if asdate:
        a = masterdate
        masterdate = dt.date(int(a[:4]),int(a[4:6]),int(a[6:8]))
    if asdatetime:
        a = masterdate
        centime = misc.grep1line('center_time',metafile)
        if not centime:
            print('error parsing center_time information from metadata.txt')
            return False
        centime = centime.split('=')[1].split('.')[0]
        masterdate = dt.datetime(int(a[:4]),int(a[4:6]),int(a[6:8]),
                        int(centime.split(':')[0]),
                        int(centime.split(':')[1]),
                        int(centime.split(':')[2]))
    return masterdate


def get_frame_master_s1ab(frame, metafile = None):
    """ Gets information if the reference epoch of given frame is S1 'A' or 'B'.
    Args:
        frame (str)
        metafile (str): if None, it will identify it on LiCSAR_public
    """
    tr = int(frame[:3])
    if not metafile:
        metafile = os.path.join(os.environ['LiCSAR_public'], str(tr), frame, 'metadata', 'metadata.txt')
    if not os.path.exists(metafile):
        print('metadata file does not exist for frame '+frame)
        return 'X'
    primepoch = grep1line('master=',metafile).split('=')[1]
    path_to_slcdir = os.path.join(os.environ['LiCSAR_procdir'], str(tr), frame, 'SLC', primepoch)
    try:
        out = os.path.basename(glob.glob(path_to_slcdir+'/S1*')[0])[2]
    except:
        print('error getting the value for frame '+frame)
        out = 'X'
    return out


def vis_aoi(aoi):
    """to visualize a polygon element ('aoi')"""
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
    plt.show()

def vis_bidtanxs(bidtanxs):
    """Visualize list of bursts (use bidtanx id, i.e. e.g. '73_IW1_1234')"""
    tovis = []
    for bid in bidtanxs:
        tovis.append(lq.get_polygon_from_bidtanx(bid))
    vis_aoi(tovis)


def vis_frame(frame):
    """Visualize frame ID"""
    ai = lq.get_bursts_in_frame(frame)
    bidtanxs = lq.sqlout2list(ai)
    vis_bidtanxs(bidtanxs)

def extract_bursts_by_track(bidtanxs, track):
    newbids = []
    for bidtanx in bidtanxs:
        if bidtanx.split('_')[0] == str(track):
            newbids.append(bidtanx)
    return newbids


def bursts2geopandas(bidtanxs, merge = False, use_s1burst = False):
    # in order to export to KML:
    # frame_gpd.to_file('~/kmls/'+frame+'.kml', driver='KML')
    # or to SHP:
    # frame_gpd.to_file('~/shps/'+frame+'.shp', driver='ESRI Shapefile')
    geometry = []
    crs = {'init': 'epsg:4326'}
    orbdir = lq.get_orbdir_from_bidtanx(bidtanxs[0])
    if merge == False:
        #if use_s1burst:
        if type(bidtanxs)==list:
            for bid in bidtanxs:
                if use_s1burst:
                    geometry.append(lq.get_s1b_geom_from_bidtanx(bid, opass = orbdir))
                else:
                    geometry.append(lq.get_polygon_from_bidtanx(bid))
            df_name = {'burstID': bidtanxs}
            aoi_gpd = gpd.GeoDataFrame(df_name, crs=crs, geometry=geometry)
        else:
            if use_s1burst:
                geometry.append(lq.get_s1b_geom_from_bidtanx(bidtanxs, opass = orbdir))
            else:
                geometry.append(lq.get_polygon_from_bidtanx(bidtanxs))
            aoi_gpd = gpd.GeoDataFrame(index=[0], crs=crs, geometry=[geometry])
    else:
        polygon = generate_frame_polygon(bidtanxs, orbdir)
        framename = generate_frame_name(bidtanxs)
        #aoi_gpd = gpd.GeoDataFrame(index=[0], crs=crs, geometry=[polygon])
        aoi_gpd = gpd.GeoDataFrame({'frameID': [framename]}, crs=crs, geometry=[polygon])
    return aoi_gpd


def frame2geopandas(frame, brute = False, use_s1burst = False, merge = False):
    if use_s1burst:
        bidtanxs=lq.get_bidtanxs_in_frame(frame)
        bidtanxs=lq.sqlout2list(bidtanxs)
        return bursts2geopandas(bidtanxs, merge = merge, use_s1burst = use_s1burst)
    if brute:
        gpan = frame2geopandas_brute(frame)
    else:
        if not lq.is_in_polygs2geom(frame):
            print('frame {} has no record in polygs2geom. Recreating'.format(frame))
            gpan = frame2geopandas_brute(frame)
        else:
            #geometry = []
            crs = {'init': 'epsg:4326'}
            wkt = lq.geom_from_polygs2geom(frame)
            geom = loads(wkt)
            #gpan['frameID'] = frame
            gpan = gpd.GeoDataFrame({'frameID': [frame]}, crs=crs, geometry=[geom])
    return gpan


def frame2geopandas_brute(frame):
    bidtanxs = lq.get_bursts_in_frame(frame)
    if not bidtanxs:
        #try it once again
        bidtanxs = lq.get_bursts_in_frame(frame)
        if not bidtanxs:
            print('the frame '+frame+' is not connected to any bursts. removing the frame')
            lq.delete_frame_only(frame)
            return None
    bidtanxs = lq.sqlout2list(bidtanxs)
    try:
        newname = generate_frame_name(bidtanxs)
    except:
        print('some problem generating frame name from the bursts of frame: '+frame)
        return None
    #get the geopandas record
    gpan = bursts2geopandas(bidtanxs, merge = True)
    if gpan.empty:
        return None
    if not newname:
        return None
    if (frame[-6:] != newname[-6:]) or (frame[3] != newname[3]):
        #print('WARNING! This frame changed its definition')
        #print('{0} ---> {1}'.format(frame,newname))
        #print('framecare_rename.sh {0} {1}'.format(frame,newname))
        #rename it in database
        lq.rename_frame(frame,newname)
        #print('now we should do: rename_frame_main(frame,newname)')
        rename_frame_main(frame,newname)
    else:
        #keep the original name if bursts did not change..
        gpan['frameID']=frame
    #outgpd = outgpd.append(gpan, ignore_index=True)
    return gpan

def rename_frame_main(framename,newname, reportcsv = '/gws/nopw/j04/nceo_geohazards_vol1/public/LiCSAR_products/frameid_changes.txt'):
    """
    this function will physically rename a frame (and move folders etc.) - oh, but it doesn't touch the frame def in database! for this, use lq.rename_frame
    """
    track = str(int(framename[0:3]))
    pubpath = os.path.join(pubdir,track,framename)
    procpath = os.path.join(procdir,track,framename)
    newpubpath = os.path.join(pubdir,track,newname)
    newprocpath = os.path.join(procdir,track,newname)
    if os.path.exists(pubpath):
        os.rename(pubpath, newpubpath)
        for fileext in ['.geo.E.tif','.geo.N.tif','.geo.hgt.tif','.geo.U.tif','-poly.txt']:
            oldfile = os.path.join(newpubpath,'metadata',framename+fileext)
            newfile = os.path.join(newpubpath,'metadata',newname+fileext)
            if os.path.exists(oldfile):
                os.rename(oldfile,newfile)
    if os.path.exists(procpath):
        os.rename(procpath, newprocpath)
    if os.path.exists(os.path.join(newprocpath,framename+'-poly.txt')):
        os.remove(os.path.join(newprocpath,framename+'-poly.txt'))
    print('frame {0} renamed to {1}'.format(framename,newname))
    if not os.path.exists(reportcsv):
        with open(reportcsv, 'w') as f:
            f.write('oldname,newname\n')
    with open(reportcsv, 'a') as f:
        f.write('{0},{1}\n'.format(framename, newname))

def get_number_of_ifgs(framename):
    pubdir = os.environ['LiCSAR_public']
    track = str(int(framename[0:3]))
    pubpath = os.path.join(pubdir,track,framename)
    if not os.path.exists(pubpath):
        return 0
    ifgspath = os.path.join(pubpath,'interferograms')
    if not os.path.exists(ifgspath):
        return 0
    filenumber = len(glob.glob1(ifgspath,'2???????_2???????'))
    return filenumber


def get_epochs(framename, return_mli_tifs = False):
    pubdir = os.environ['LiCSAR_public']
    track = str(int(framename[0:3]))
    pubpath = os.path.join(pubdir,track,framename)
    if not os.path.exists(pubpath):
        return False
    epochspath = os.path.join(pubpath,'epochs')
    if not os.path.exists(epochspath):
        return False
    if return_mli_tifs:
        #this will return tif file paths
        return glob.glob(epochspath + "/**/*.geo.mli.tif", recursive = True)
    else:
        epochslist = glob.glob1(epochspath,'2???????')
        return epochslist


def get_ifg_list_pubdir(framename):
    pubdir = os.environ['LiCSAR_public']
    track = str(int(framename[0:3]))
    pubpath = os.path.join(pubdir,track,framename)
    if not os.path.exists(pubpath):
        return 0
    ifgspath = os.path.join(pubpath,'interferograms')
    if not os.path.exists(ifgspath):
        return 0
    ifglist = glob.glob1(ifgspath,'2???????_2???????')
    return ifglist


def get_epochs_from_ifg_list_pubdir(framename):
    ifglist = get_ifg_list_pubdir(framename)
    epochs = set()
    if ifglist == 0:
        return 0
    else:
        for ifg in ifglist:
            epochs.add(ifg.split('_')[0])
            epochs.add(ifg.split('_')[1])
    return list(epochs)


def export_frames_to_licsar_csv(framesgpd, outcsv = '/gws/nopw/j04/nceo_geohazards_vol1/public/shared/frames/frames.csv', store_zero = False):
    #print('now we would export the frame to outcsv, including wkb')
    # this will update the csv, not rewrite it..
    if not os.path.exists(outcsv):
        with open(outcsv,'w') as f:
            f.write('the_geom,frame,files,download,direction\n')
    with open(outcsv,'a') as f:
        #extract needed information
        for index, row in framesgpd.iterrows():
            framename = row['frameID']
            track = int(framename[0:3])
            download = "<a href='http://gws-access.ceda.ac.uk/public/nceo_geohazards/LiCSAR_products/{0}/{1}/' target='_blank'>Link<a>".format(str(track),framename)
            orbdir = framename[3]
            if orbdir == 'A':
                direction = 'Ascending'
            elif orbdir == 'D':
                direction = 'Descending'
            else:
                print('wrong framename! Aborting')
                return False
            #get number of files
            files = get_number_of_ifgs(framename)
            #get geometrywkb
            geom = row['geometry'].wkb_hex
            #we do not want to export it in case of no files in public...:
            if ((not store_zero) and files > 0) or store_zero:
                #if frame already in csv file, remove its line and update by new files no.
                if misc.grep1line(framename, outcsv):
                    misc.sed_rmlinematch(framename, outcsv)
                f.write('{0},{1},{2},{3},{4}\n'.format(str(geom), framename, str(files), download, direction))


def store_frame_geometry(framesgpd):
    fileio_error = False
    for index, row in framesgpd.iterrows():
        framename = row['frameID']
        track = int(framename[0:3])
        geom = row['geometry'].wkt
        #update the xy file:
        pubfile = os.path.join(pubdir,str(track),framename,'metadata',framename+'-poly.txt')
        procfile = os.path.join(procdir,str(track),framename,framename+'-poly.txt')
        procframefile = os.path.join(procdir,str(track),framename,'frame.xy')
        xy = row['geometry'].exterior.coords.xy
        x = xy[0]
        y = xy[1]
        for fileout in [pubfile, procfile, procframefile]:
            if os.path.exists(fileout):
                os.remove(fileout)
            try:
                with open(fileout,'w') as f:
                    for i in range(len(x)-1):
                        f.write(str(x[i])+' '+str(y[i])+'\n')
            except:
                fileio_error = True
                #print('warning, {0} could not have been generated'.format(fileout))
        #update the database GIS table
        res = lq.store_frame_geometry(framename, geom)
    return res

def export_geopandas_to_kml(gpan, outfile):
    gpan.to_file(outfile, driver='KML', NameField='frameID')

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

def generate_frame_polygon(bidtanxs, orbdir = None):
    if not orbdir:
        orbdir = lq.get_orbdir_from_bidtanx(bidtanxs[0])
    try:
        burstgpd = bursts2geopandas(bidtanxs)
    except:
        print('some error during bursts2geopandas, maybe mysql problem')
        return None
    #unite bursts, but this will keep errors:
    framegpd = burstgpd.unary_union
    #corrections based on:
    # https://gis.stackexchange.com/questions/277334/shapely-polygon-union-results-in-strange-artifacts-of-tiny-non-overlapping-area
    eps = 0.025
    tolsim = eps
    from shapely.geometry import JOIN_STYLE
    framegpd = framegpd.buffer(eps, 1, join_style=JOIN_STYLE.mitre).buffer(-eps, 1, join_style=JOIN_STYLE.mitre)
    framegpd = framegpd.simplify(tolerance=tolsim)
    #maximal number of points should be 13!
    while len(framegpd.exterior.coords[:])>13:
        tolsim = tolsim+0.001
        framegpd = framegpd.simplify(tolerance=tolsim)
    return Polygon(framegpd.exterior)

def generate_frame_polygon_old(bidtanxs, orbdir):
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
    orbdir = lq.get_orbdir_from_bidtanx(bidtanxs[0])
    polyhon = generate_frame_polygon(bidtanxs, orbdir)
    if not polyhon:
        print('some error generating frame polygon - mysql access error?')
        return None
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

def generate_new_frame(bidtanxs,testonly = True, hicode = None):
    """
    hicode... use 'H' for high resolution (1/5 multilook), or 'M' for medium, i.e. 56 m
    """
    #and now i can generate the new frame:
    track = bidtanxs[0].split('_')[0]
    orbdir = lq.get_orbdir_from_bidtanx(bidtanxs[0])
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
    if hicode:
        polyid_name = polyid_name[:9]+hicode+polyid_name[10:]
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
    try:
        polyid = int(lastpolyid[0][0])+1
    except:
        print('seems like first frame to ingest - ok')
        polyid = 1
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
        print('including to polyg2gis table')
        gpan = frame2geopandas_brute(polyid_name)
        rc = store_frame_geometry(gpan)
        if rc != 1:
            print('ERROR STORING TO polyg2gis TABLE!!!')
        #else:
        #
        #actually the licsar csv should not contain this..
        #rc = export_frames_to_licsar_csv(gpan)
        print('generated new frame '+polyid_name)
        print('you may do following now: ')
        print('licsar_initiate_new_frame.sh '+polyid_name)
        if hicode:
            print('(but remember adding parameter -'+hicode+')')
        #delete_frame_commands(frame)
    return polyid_name


def generate_frame_from_bursts_kml(inputkml):
    bursts = load_bursts_from_kml(inputkml)
    generate_new_frame(bursts, testonly=False)


def load_bursts_from_txt(intxt):
    f = open(intxt,'r')
    contents = f.readlines()
    bursts = []
    for burst in contents:
        bursts.append(burst.split('\n')[0])
    f.close()
    return bursts


def load_xy(intxt, onlyminmax = True):
    f = open(intxt,'r')
    contents = f.readlines()
    f.close()
    lon = []
    lat = []
    for line in contents:
        lon.append(float(line.split(' ')[0]))
        lat.append(float(line.split(' ')[1]))
    if onlyminmax:
        out = [min(lon), max(lon)], [min(lat), max(lat)]
    else:
        out = [lon, lat]
    return out


def load_bursts_from_kml(inputkml):
    #inputkml = '/gws/nopw/j04/nceo_geohazards_vol2/LiCS/temp/insar_temp/frames_redef/kmls/test_146a.kml'
    newbursts = gpd.read_file(inputkml, driver='KML')
    newbursts = newbursts[newbursts.columns[0]].tolist()
    return newbursts



def export_bidtanxs_to_kml(bidtanxs, outpath = '/gws/nopw/j04/nceo_geohazards_vol1/public/shared/test', projname = 'track', merge = False):
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

def export_frame_to_kml(frame, outpath = '/gws/nopw/j04/nceo_geohazards_vol2/LiCS/temp/insar_temp/frames_redef/kmls', merge=False):
    """ Exports the frame polygon to kml. Currently only the uglier version (coarse burst polygons)
    """
    if not os.path.exists(outpath):
        os.mkdir(outpath)
    ai = lq.get_bursts_in_frame(frame)
    bidtanxs = lq.sqlout2list(ai)
    export_bidtanxs_to_kml(bidtanxs, outpath, projname = frame, merge=merge)


def export_all_frames_to_framecsv(outcsv = '/gws/nopw/j04/nceo_geohazards_vol1/public/shared/frames/frames.csv', store_zero = False):
    asc_gpd, desc_gpd = get_all_frames()
    rc = export_frames_to_licsar_csv(asc_gpd, outcsv, store_zero)
    rc = export_frames_to_licsar_csv(desc_gpd, outcsv, store_zero)
    return rc


def get_satids_of_burst(burstid, expected = ['S1A','S1B']):
    a = lq.get_filenames_from_burst(burstid)
    a = lq.sqlout2list(a)
    b = pd.DataFrame(expected)
    b['count'] = 0
    for x in a:
        satid=x[:3]
        if satid in expected:
            i=b[b[0] == satid]['count'].index[0]
            b.loc[i,'count'] =b.loc[i,'count']+1
    b = b.rename(columns={0: "sat_id"})
    return b


def get_all_frames():
    asc_gpd = gpd.geodataframe.GeoDataFrame()
    desc_gpd = gpd.geodataframe.GeoDataFrame()
    for i in range(1,175+1):
        print('preparing frames from track {}'.format(i))
        #descending:
        frames = lq.get_frames_in_orbit(i, 'D')
        frames = lq.sqlout2list(frames)
        for frame in frames:
            a = frame2geopandas(frame)
            if type(a) != type(None):
                desc_gpd = desc_gpd.append(a)
        #ascending
        frames = lq.get_frames_in_orbit(i, 'A')
        frames = lq.sqlout2list(frames)
        for frame in frames:
            a = frame2geopandas(frame)
            if type(a) != type(None):
                asc_gpd = asc_gpd.append(a)
    return asc_gpd, desc_gpd


def manual_check_master_files(frame, master):
    if type(master) == type('str'):
        masterdate=dt.datetime.strptime(master,'%Y%m%d')
    else:
        masterdate=master
    filelist = lq.get_frame_files_date(frame,masterdate)
    print(len(filelist))
    for filee in filelist:
        fid=filee[1]
        brsts = lq.get_bursts_in_file(fid)
        brsts = lq.sqlout2list(brsts)
        b=bursts2geopandas(brsts)
        print(filee[2])
        b.plot()
        plt.show()
    print('to fix the bad ones:')
    print('fullpath = ....')
    print('lq.delete_file_from_db(fullpath)')
    print("os.system('arch2DB.py -f {} >/dev/null 2>/dev/null'.format(fullpath))")


def export_all_frames_to_kmls(kmldirpath = '/gws/nopw/j04/nceo_geohazards_vol1/public/shared/frames/'):
    asc_gpd, desc_gpd = get_all_frames()
    if os.path.exists(os.path.join(kmldirpath,'ascending.kml')):
        os.remove(os.path.join(kmldirpath,'ascending.kml'))
    if os.path.exists(os.path.join(kmldirpath,'descending.kml')):
        os.remove(os.path.join(kmldirpath,'descending.kml'))
    
    export_geopandas_to_kml(asc_gpd, os.path.join(kmldirpath,'ascending.kml'))
    export_geopandas_to_kml(desc_gpd, os.path.join(kmldirpath,'descending.kml'))


def delete_bursts(bidtanxs, test = True):
     if test:
         print('WARNING, this will delete all bursts in the list FOREVER')
         print('if any frame still uses the burst ids, it will cancel their deletion')
         return
     else:
         for bidtanx in bidtanxs:
             print('deleting burst '+bidtanx)
             try:
                 rc = lq.delete_burst_from_db(bidtanx)
             except:
                 print('some error occurred - cancelling')
                 return


def delete_frame(frame):
    #print('cannot use this anymore, in CentOS7 - please contact admin to delete frame '+frame)
    #return False
    polyid = lq.get_frame_polyid(frame)[0][0]
    if not polyid:
        print('error - is it correct frame??')
        return
    #cmd_mysql='mysql -h ..... licsar_batch ' # see lics_mysql.sh
    os.system('setFrameInactive.py {0}'.format(frame))
    track=str(int(frame[0:3]))
    os.system('rm -rf $LiCSAR_procdir/{0}/{1} $LiCSAR_public/{0}/{1}'.format(track,frame))
    os.system("sed -i '/{}/d' /gws/nopw/j04/nceo_geohazards_vol1/public/shared/frames/frames.csv".format(frame))
    os.system("sed -i '/{}/d' /gws/nopw/j04/nceo_geohazards_vol1/public/LiCSAR_products/EQ/eqframes.csv".format(frame))
    #polyid = lq.get_frame_polyid(frame)[0][0]
    sql = "delete from licsar_batch.polygs2master where polyid={};".format(polyid)
    rc = lq.do_query(sql, 1)
    #os.system(cmd_mysql+' -e "{}"'.format(sql))
    sql = "delete from polygs2bursts where polyid={};".format(polyid)
    rc = lq.do_query(sql, 1)
    try:
        sql = "delete from esd where polyid={};".format(polyid)
        rc = lq.do_query(sql, 1)
    except:
        print('error in deleting from esd table')
    #os.system(cmd_mysql+' -e "{}"'.format(sql))
    frame_workaround = frame.replace('A','Y')
    frame_workaround = frame_workaround.replace('D','Y')
    sql = "update polygs set polyid_name='{0}' where polyid_name='{1}';".format(frame_workaround, frame)
    rc = lq.do_query(sql, 1)
    #os.system(cmd_mysql+' -e "{}"'.format(sql))
    sql = "delete from polygs where polygs.polyid_name='{}';".format(frame_workaround)
    rc = lq.do_query(sql, 1)
    #rc = os.system(cmd_mysql+' -e "{}" 2>/dev/null'.format(sql))
    #if rc != 0:
    #    print('WARNING: the frame was only partially removed. But it should not appear in processing')
    print('the frame {} has been removed, associated files purged'.format(frame))


def delete_frame_commands(frame):
    print('setFrameInactive.py {0}'.format(frame))
    track=str(int(frame[0:3]))
    print('rm -rf $LiCSAR_procdir/{0}/{1} $LiCSAR_public/{0}/{1}'.format(track,frame))
    print("sed -i '/{}/d' /gws/nopw/j04/nceo_geohazards_vol1/public/shared/frames/frames.csv".format(frame))
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

def remove_bad_bursts(frame, badbursts, testonly = True):
    #badbursts should be a list of bursts existing in the frame that should be removed
    #e.g. the list of missing bursts during the licsar_init_frame.sh script..
    #in testonly - it will only give text output, rather than really do something..
    bursts = lq.get_bidtanxs_in_frame(frame)
    bursts = lq.sqlout2list(bursts)
    for b in badbursts:
        bursts.remove(b)
    generate_new_frame(bursts, testonly)
    print('to remove the old frame, you should do:')
    print("fc.delete_frame('{}')".format(frame))
    #delete_frame_commands(frame)


def add_more_bursts(frame, extrabursts, testonly = True):
    #extrabursts should be a list of bursts not existing in the frame, to be added there
    #in testonly - it will only give text output, rather than really do something..
    bursts = lq.get_bidtanxs_in_frame(frame)
    bursts = lq.sqlout2list(bursts)
    newbursts = bursts+extrabursts
    #remove duplicities
    newbursts = list(set(newbursts))
    newbursts.sort()
    generate_new_frame(newbursts, testonly)
    print('to remove the old frame, you should do:')
    print("fc.delete_frame('{}')".format(frame))


