#!/usr/bin/env python

import os, shutil
import LiCSquery as lq
import datetime as dt
from datetime import datetime, timedelta
from libcomcat.search import search,count,get_event_by_id
import pandas as pd
from LiCSAR_lib.s1data import get_images_for_frame
import LiCSAR_lib.framecare as fc
import numpy as np
import LiCSAR_lib.LiCSAR_misc as misc
import framecare as fc
import time

public_path = os.environ['LiCSAR_public']
procdir_path = os.environ['LiCSAR_procdir']
web_path = 'https://gws-access.ceda.ac.uk/public/nceo_geohazards/LiCSAR_products'
web_path_maps = 'https://comet.nerc.ac.uk/earthquakes'

eqcsvfile = "/home/home02/earmla/pokuseq.csv"

#you may want to change these parameters:
max_days = 35
minmag = 5.5
#here will be exceptions (i.e. eqs that MUST be processed):
#minmag = 4.6
exceptions = [] #'us60008e8e', 'us6000a50y']


#table based on John Elliott's know-how..
eq_data = {'magnitude': [5.5,5.6,5.7,5.8,5.9,6, \
                              6.1,6.2,6.3,6.4,6.5,6.6,6.7,6.8,6.9,7, \
                              7.1,7.2,7.3,7.4,7.5,7.6,7.7,7.8,7.9,8, \
                              8.1,8.2,8.3,8.4,8.5],
          'distance': [20,20,21,25,30,36, \
                      42,50,60,71,84,100,119,141,167,199, \
                      236,280,333,396,471,559,664,790,938,1115, \
                      1325, 1575, 1872, 2225, 2500],
          'depth': [10,11,12,14,16,18, \
                   21,24,28,32,36,42,48,55,63,73, \
                   83,96,110,127,146,168,193,222,250,250, \
                   250, 250, 250, 250, 250]}
eq_limits = pd.DataFrame(eq_data, columns = {'magnitude', 'distance', 'depth'})



def get_range_from_magnitude(M, depth, unit = 'km'):
    #M = 5.5
    #depth = 9 #km
    M = round(M,1)
    if M > 8.5 and depth <= 250:
        distance = 2500
    else:
        distance = eq_limits.query('magnitude == {0} and depth >= {1}'.format(M,depth))['distance']
        if len(distance)>0:
            distance = int(distance.to_string().split()[1])
        else:
            distance = None
            #print('The earthquake parameters do not fit with the limit table - it will not be processed')
    #rad or km?
    if unit == 'rad' and distance is not None:
        distance = (360.0 / 40007.86) * distance
    return distance

def create_eq_csv(csvfile = '/gws/nopw/j04/nceo_geohazards_vol1/public/LiCSAR_products/EQ/eqs.csv'):
    #get all eqs
    query = "select * from eq;"
    eq_df = lq.do_pd_query(query)
    eq_df['link'] = "<a href='{}/".format(web_path_maps) + eq_df['USGS_ID'] + ".html' target='_blank'>Link</a>"
    dbcols = ['USGS_ID','magnitude','depth','time','lat','lon', 'link', 'location']
    eq_df[dbcols].to_csv(csvfile, sep = ';', index=False)
    return True

def add_new_event(event, updatecsv = True):
    #this will add the event to the database of eqs and the eqs.csv
    if type(event)==str:
        #assuming eventid is parsed only
        event =  get_event_by_id(event)
    
    #let's update database
    eqid = import_to_licsinfo_eq(event)
    if not eqid:
        #just update it then
        lq.do_query("update eq set location='{0}' where USGS_ID='{1}'".format(event.location.replace("'"," "), event.id), 1)
    if updatecsv:
        update_eq_csv(event.id, eqcsvfile)
    return True


def update_eq_csv(eventid, csvfile = '/gws/nopw/j04/nceo_geohazards_vol1/public/LiCSAR_products/EQ/eqs.csv', delete=False):
    if os.path.exists(csvfile+'.lock'):
        print('file lock detected, waiting 10 minutes - should be enough for ARC to finish')
        time.sleep(10*60)
    if os.path.exists(csvfile+'.lock'):
        print('probably error during copy by ARC_responder. removing lock file')
        os.remove(csvfile+'.lock')
    if not delete:
        query = "select * from eq where USGS_ID='{}';".format(eventid)
        eq_df = lq.do_pd_query(query)
        eq_df['link'] = "<a href='{}/".format(web_path_maps) + eq_df['USGS_ID'] + ".html' target='_blank'>Link</a>"
        dbcols = ['USGS_ID','magnitude','depth','time','lat','lon', 'link', 'location']
        eq_df[dbcols].to_csv(csvfile, mode='a', sep = ';', index = False, header=False)
        return True
    else:
        rc = os.system("sed -i '/{0}/d' {1}".format(eventid, csvfile))
        if rc == 0:
            return True
        else:
            print('error during removal of the eventid '+eventid)
            return False


def regenerate_eq2frames_csv(csvfile = '/gws/nopw/j04/nceo_geohazards_vol1/public/LiCSAR_products/EQ/eqframes.csv'):
    query = "select p.polyid_name as frame, aswkb(pg.geom) as the_geom, eq.USGS_ID as usgsid, eq.location, e2f.post_acq as next_possible, e2f.next_acq as next_expected  \
        from eq2frame e2f inner join eq on e2f.eqid=eq.eqid inner join polygs2gis pg  \
        on pg.polyid=e2f.fid inner join polygs p on p.polyid=e2f.fid"
    e2f = lq.do_pd_query(query)
    if len(e2f) < 1:
        print('error')
        return False
    e2f['track'] = e2f.frame.apply(lambda x:str(int(x[:3])))
    e2f['download'] = e2f['track']*0
    for i, f in e2f.iterrows():
        #check if we have geometry here:
        if len(f.the_geom) < 5:
            try:
                row = fc.frame2geopandas(f['frame']).iloc[0]
                geometry = row['geometry'].wkb_hex
                e2f.at[i, 'the_geom'] = geometry
            except:
                print('ERROR: probably frame {0} is not in geom database'.format(f['frame']))
                continue
        else:
            #fix the strange format to ordinary wkb_hex
            if str(e2f.at[i, 'the_geom'])[0] == 'b':
                e2f.at[i, 'the_geom'] = e2f.at[i, 'the_geom'].hex().upper()
        downlink = "<a href='{0}/{1}/{2}/' target='_blank'>Link</a>".format(web_path, f['track'], f['frame'])
        #if we do not have this combination of frame and event id, include it
        e2f.at[i, 'download'] = downlink
    #e2f['download'] = "<a href='http://gws-access.ceda.ac.uk/public/nceo_geohazards/LiCSAR_products/{0}/{1} target='_blank'>Link<a>".format(track, framename)
    #dbcols = ['the_geom','frame','usgsid','download', 'location','next_possible', 'next_expected']
    dbcols = ['the_geom','frame','usgsid','download', 'next_possible', 'next_expected']
    e2f[dbcols].query('download != ""').to_csv(csvfile, mode='w', sep = ';', index = False, header=False)
    return True


def reset_frame(frame, eventid='us7000ckx5'):
    '''
    this will set the frame to coifg status 0 and frame status 2, i.e. copied2arc, checking for coseismic ifg
    '''
    try:
        fid=lq.get_frame_polyid(frame)[0][0]
    except:
        print('error getting frame ID')
        return False
    try:
        eqid = lq.get_eqid(eventid)
    except:
        print('error getting eq ID')
        return False
    lq.update('eq2frame','coifg_status','0','fid={} and eqid={}'.format(fid, eqid))
    lq.update('eq2frame','frame_status','2','fid={} and eqid={}'.format(fid, eqid))
    lq.update('eq2frame','postifg_status','0','fid={} and eqid={}'.format(fid, eqid))
    return True


def update_eq2frames_csv(eventid, csvfile = '/gws/nopw/j04/nceo_geohazards_vol1/public/LiCSAR_products/EQ/eqframes.csv', delete = False):
    """
       This csv will be loaded to the eq responder map
       WARNING - would read all event frames from database and only add them to the csv file if it is not there yet
    """
    dbcols = ['the_geom','frame','usgsid','download', 'next_possible', 'next_expected']
    if os.path.exists(csvfile+'.lock'):
        print('file lock detected, waiting 10 minutes - should be enough for ARC to finish')
        time.sleep(10*60)
    if os.path.exists(csvfile+'.lock'):
        print('probably error during copy by ARC_responder. removing lock file')
        os.remove(csvfile+'.lock')
    if delete:
        print('deleting of event '+eventid+' from eqframes.csv')
        rc = os.system("sed -i '/{0}/d' {1}".format(eventid, csvfile))
        if rc == 0:
            return True
        else:
            print('error during removal of the eventid '+eventid+' from eqframes.csv')
            return False
    #fid = lq.get_frame_polyid(framename)
    #try:
    #    fid = lq.sqlout2list(fid)[0]
    #except:
    #    print('error - the frame {} does not exist in licsinfo database'.format(framename))
    #query = "select pg.geom from eq2frame e2f inner join eq on e2f.eqid=eq.eqid inner join polygs2gis pg on pg.polyid=e2f.fid where eq.USGS_ID='{0}' and e2f.fid={1};".format(eventid, str(fid))
    #query = "select pg.geom,e2f.fid from eq2frame e2f inner join eq on e2f.eqid=eq.eqid inner join polygs2gis pg on pg.polyid=e2f.fid where eq.USGS_ID='{0}'".format(eventid)
    #query = "select p.polyid_name as frame, aswkb(pg.geom) as the_geom, eq.USGS_ID as usgsid, eq.location, e2f.post_acq as next_possible, e2f.next_acq as next_expected  \
    query = "select p.polyid_name as frame, aswkb(pg.geom) as the_geom, eq.USGS_ID as usgsid, eq.location, e2f.post_acq as next_possible, e2f.next_acq as next_expected  \
        from eq2frame e2f inner join eq on e2f.eqid=eq.eqid inner join polygs2gis pg  \
        on pg.polyid=e2f.fid inner join polygs p on p.polyid=e2f.fid where eq.USGS_ID='{0}'".format(eventid)
        #this is the original query... it would not ingest frames that fail in the beginning (copy from CEDA or initialisation). now keeping all identified frames..
        #on pg.polyid=e2f.fid inner join polygs p on p.polyid=e2f.fid where eq.USGS_ID='{0}' and e2f.frame_status > 0 and e2f.frame_status < 50".format(eventid)
    e2f = lq.do_pd_query(query)
    if len(e2f) < 1:
        print('error - nothing found in eq2frames for this event {}'.format(eventid))
        return False
    e2f['track'] = e2f.frame.apply(lambda x:str(int(x[:3])))
    e2f['download'] = e2f['track']*0
    for i, f in e2f.iterrows():
        #check if we have geometry here:
        if len(f.the_geom) < 5:
            try:
                row = fc.frame2geopandas(f['frame']).iloc[0]
                geometry = row['geometry'].wkb_hex
                e2f.at[i, 'the_geom'] = geometry
            except:
                print('ERROR: probably frame {0} is not in geom database'.format(f['frame']))
                return False
        else:
            #fix the strange format to ordinary wkb_hex
            if str(e2f.at[i, 'the_geom'])[0] == 'b':
                e2f.at[i, 'the_geom'] = e2f.at[i, 'the_geom'].hex().upper()
        downlink = "<a href='{0}/{1}/{2}/' target='_blank'>Link</a>".format(web_path, f['track'], f['frame'])
        #f.download = downlink
        e2f.at[i, 'download'] = downlink
        strincsv = e2f[dbcols].loc[i].to_csv(sep = ';', index = False, header=False).replace('\n',';').replace('"','')[:-1]
        #minimal string to see if there is related line (to be deleted then)
        minstrincsv = f.frame+';'+f.usgsid
        greppedstr = misc.grep1line(minstrincsv, csvfile)
        if greppedstr:
            if greppedstr == strincsv:
                #so this line exists i csvfile, so skipping it..
                e2f = e2f.drop(index=i)
            else:
                #we will need to update this line..
                rc = os.system("sed -i '/{0}/d' {1}".format(minstrincsv, csvfile))
        #if we do not have this combination of frame and event id, include it
        #if not misc.grep1line(f['frame']+';'+eventid, csvfile):
        #    e2f.at[i, 'download'] = downlink
        #else:
        #    e2f.at[i, 'download'] = ''
    #e2f['download'] = "<a href='http://gws-access.ceda.ac.uk/public/nceo_geohazards/LiCSAR_products/{0}/{1} target='_blank'>Link<a>".format(track, framename)
    #dbcols = ['the_geom','frame','usgsid','download', 'location','next_possible', 'next_expected']
    e2f[dbcols].to_csv(csvfile, mode='a', sep = ';', index = False, header=False)
    #e2f[dbcols].query('download != ""').to_csv(csvfile, mode='a', sep = ';', index = False, header=False)
    return True


def list_coseismic_ifgs(frame, toi):
    #this would list all ifgs in public directory that are coseismic (if there are such)
    #toi is datetime
    doi = toi.date()
    track = str(int(frame[0:3]))
    global public_path
    products_path = os.path.join(public_path, track, frame, 'interferograms')
    if not os.path.exists(products_path):
        print('error - frame has no interferograms')
        return False
    doi_str = int(doi.strftime('%Y%m%d'))
    publicifgs = os.listdir(products_path)
    if not publicifgs:
        print('error - frame has no interferograms')
        return []
    publicifgs = [i for i in publicifgs if '_' in i]
    selected_ifgs = []
    for pifg in publicifgs:
        is_coseismic = False
        mas = int(pifg.split('_')[0])
        slv = int(pifg.split('_')[1])
        #this is to generate kmz files only for coseismic ifgs
        if (doi_str > mas) and (doi_str < slv):
            is_coseismic = True
        if (doi_str == mas):
            date = datetime.strptime(str(mas),'%Y%m%d').date()
            filelist = lq.get_frame_files_date(frame, date)
            tof = lq.get_time_of_file(filelist[0][1])
            if tof < toi:
                is_coseismic = True
        if (doi_str == slv):
            date = datetime.strptime(str(slv),'%Y%m%d').date()
            filelist = lq.get_frame_files_date(frame, date)
            if filelist:
                tof = lq.get_time_of_file(filelist[0][1])
                if tof > toi:
                    is_coseismic = True
            else:
                is_coseismic = False
        if is_coseismic:
            selected_ifgs.append(pifg)
    return selected_ifgs


def create_kmls(frame, toi, onlycoseismic = False, overwrite = False):
    #toi is Time Of Interest - and should be as datetime
    doi = toi.date()
    track = str(int(frame[0:3]))
    global public_path
    products_path = os.path.join(public_path, track, frame, 'interferograms')
    frameproc_path = os.path.join(procdir_path, track, frame)
    if not os.path.exists(products_path):
        try:
            os.mkdir(products_path)
        except:
            print('error - this directory could not have been created: '+products_path)
            return None
    doi_str = int(doi.strftime('%Y%m%d'))
    publicifgs = os.listdir(products_path)
    publicifgs = [i for i in publicifgs if '_' in i]
    if publicifgs is None:
        print('there are no final georeferenced products generated')
        return None
    ## should do check for valid folder name here!
    selected_ifgs = []
    for pifg in publicifgs:
        is_coseismic = False
        is_postseismic = False
        is_preseismic = False
        mas = int(pifg.split('_')[0])
        slv = int(pifg.split('_')[1])
        #this is to generate kmz files only for coseismic ifgs
        if (doi_str > mas) and (doi_str < slv):
            is_coseismic = True
        if (doi_str < mas) and (doi_str < slv):
            is_postseismic = True
        if (doi_str > mas) and (doi_str > slv):
            is_preseismic = True
        if (doi_str == mas):
            date = datetime.strptime(str(mas),'%Y%m%d').date()
            filelist = lq.get_frame_files_date(frame, date)
            tof = lq.get_time_of_file(filelist[0][1])
            if tof < toi:
                is_coseismic = True
            if tof > toi:
                is_postseismic = True
        if (doi_str == slv):
            date = datetime.strptime(str(slv),'%Y%m%d').date()
            #print('debug')
            #print(date)
            #print('file list is:')
            filelist = lq.get_frame_files_date(frame, date)
            if filelist:
                tof = lq.get_time_of_file(filelist[0][1])
                if tof > toi:
                    is_coseismic = True
                if tof < toi:
                    is_preseismic = True
            else:
                is_coseismic = False
                is_preseismic = False
        if onlycoseismic:
            cond = is_coseismic
        else:
            cond = (is_coseismic or is_postseismic)
        if cond:
            if not overwrite:
                if not os.path.exists(os.path.join(products_path,pifg,frame+'_'+pifg+'.kmz')):
                    os.system('create_kmz.sh {0} {1}'.format(os.path.join(products_path,pifg), frame))
                else:
                    if dt.datetime.fromtimestamp(os.path.getctime(os.path.join(products_path,pifg,frame+'_'+pifg+'.kmz'))) <= pd.to_datetime('2020-11-01'):
                        print('detected old version of the kmz - overwriting it')
                        os.system('create_kmz.sh {0} {1}'.format(os.path.join(products_path,pifg), frame))
            else:
                os.system('create_kmz.sh {0} {1}'.format(os.path.join(products_path,pifg), frame))
            if not os.path.exists(os.path.join(products_path,pifg,frame+'_'+pifg+'.kmz')):
                print('some error happened during KMZ generation of '+pifg)
            else:
                selected_ifgs.append(pifg)
    return selected_ifgs

def get_earliest_expected_dt(frame, eventtime, metafile = None):
    if metafile:
        masterdate = fc.get_master(frame, asdatetime = True, metafile = metafile)
    else:
        masterdate = fc.get_master(frame, asdatetime = True)
    if not masterdate:
        print('error getting masterdate')
        return False
    daysdiff = (eventtime - masterdate).days
    noepochs = int(np.floor(daysdiff/6))
    lastepoch = masterdate+dt.timedelta(days=noepochs*6)
    expected_dt = lastepoch+dt.timedelta(days=6)
    return expected_dt


"""
john's table about post_acq last date:
<6.0    -->   3         (typically 18-36 post eq days)
6.0-6.5     -->   5   (typically 30-60 post eq days)
6.5-7.0     -->   7   (typically 42-84 post eq days)
7.0-7.5     -->   9   (typically 54-108 post eq days)
7.5-8.0     -->   11  (typically 66-132 post eq days, i.e. 2-4 months)
8+             -->   13  (typically 108-156 post eq days, i.e. 3-6 months)

a workaround for dummy eq - if M0 == last date will be after... X acquisitions..
"""
def get_max_post_days(magnitude, daysdiff):
    x = 100 #should give some two years...
    if magnitude == 0:
        noimag = x
    elif magnitude < 6:
        noimag = 3
    elif magnitude < 6.5:
        noimag = 5
    elif magnitude < 7:
        noimag = 7
    elif magnitude < 7.5:
        noimag = 9
    elif magnitude < 8:
        noimag = 5
    else:
        noimag = 13
    days = noimag*daysdiff
    return days


def get_next_expected_datetime(frame, eventtime, revisit_days = 6):
    #returns next expected image date, earliest possible date and expected day difference (revisit time)
    track = str(int(frame[0:3]))
    lastmonth = eventtime-timedelta(days=30)
    eqimages = get_images_for_frame(frame, startdate = lastmonth.date(), enddate = eventtime.date(), sensType = 'IW')
    eqimgdates = set()
    for eqi in eqimages:
        imgdate = eqi.split('T')[0].split('_')[-1]
        #imgtime = eqi.split('T')[1].split('_')[0]
        #imgdatetime = imgdate+'T'+imgtime
        #imgtime = datetime.strptime(imgtime,'%H%M%S').time()
        imgdate = datetime.strptime(imgdate,'%Y%m%d')
        eqimgdates.add(imgdate)
    eqimgdates = sorted(eqimgdates) # list now
    
    #lastimage = eqimgdates[-1]
    expectedtime = eqimages[-1].split('T')[1].split('_')[0]
    expectedtime = datetime.strptime(expectedtime,'%H%M%S').time()
    nextposimage = eqimgdates[-1]+timedelta(days=revisit_days)
    nextposimage = nextposimage.replace(hour=expectedtime.hour, minute=expectedtime.minute, second=expectedtime.second)
    
    lag = 0
    while nextposimage < eventtime:
        lag = lag+1
        nextposimage = nextposimage+timedelta(days=revisit_days)
    if lag>0:
        print('warning, there was no acquisition from {} pre-event pass(es)'.format(str(lag)))
    
    imgdiffs = []
    for i in range(len(eqimgdates)-1):
        imgdiffs.append((eqimgdates[i+1] - eqimgdates[i]).days)
    bestcasediff = min(imgdiffs)
    worstcasediff = max(imgdiffs)
    if not worstcasediff == bestcasediff:
        print('Temporal sampling is irregular (expecting worst case), i.e. '+str(worstcasediff)+' days diff')
    imgdiff = worstcasediff
    
    nextexpectedimage = eqimgdates[-1] + timedelta(days=imgdiff)
    nextexpectedimage = nextexpectedimage.replace(hour=expectedtime.hour, minute=expectedtime.minute, second=expectedtime.second)
    while nextexpectedimage < eventtime:
        nextexpectedimage = nextexpectedimage+timedelta(days=revisit_days)
    
    returnlist = [nextexpectedimage, nextposimage, imgdiff]
    return returnlist

def get_next_expected_images(frames, eventtime):
    print('checking first expected post-seismic acquisition time for covering frames:')
    for frame in frames:
        frame = frame[0]
        print(frame+': ')
        [nextexp, nextpos, daydiff] = get_next_expected_datetime(frame, eventtime)
        print('Expected: '+str(nextexp))
        if nextexp != nextpos:
            print('Warning, this area is not observed at the highest frequency')
            print('Next possible flight: '+str(nextpos))

def get_eq_events(minmag = 5.5, max_days = 30):
    try:
        out = search(starttime=datetime.now()-timedelta(days=max_days),
                           endtime=datetime.now(),
                           minmagnitude=minmag)
    except:
        out = False
    return out

def get_event_details(eventcode):
    eventlist = get_eq_events()
    for event in eventlist:
        if event.id == eventcode:
            selevent = event
    if selevent:
        return selevent
    else:
        return None

def get_frames_in_event(event,radius = 9999):
    if radius == 9999:
        radius = get_range_from_magnitude(event.magnitude, event.depth, 'rad')
    minlon = event.longitude - radius
    maxlat = event.latitude + radius
    maxlon = event.longitude + radius
    minlat = event.latitude - radius
    frames = lq.get_frames_in_polygon(minlon,maxlon,minlat,maxlat)
    return frames


def import_to_licsinfo_eq(event, active = True):
    if lq.get_eqid(event.id):
        print('eq already in database')#', cancelling')
        return False
    else:
        eqid = lq.insert_new_eq(event, active)
        return eqid


def import_to_licsinfo_eq2frame(eqid, event, frame, postacq = True, active = True):
    '''
    checks and inserts frame of a particular event using insert_new_eq2framme
    '''
    fid = lq.get_frame_polyid(frame)
    try:
        fid = lq.sqlout2list(fid)[0]
    except:
        print('the frame does not exist in licsinfo!')
        return False
    rc = False
    #in case this is the dummy event (i.e. should use the current time)
    #if eqid == 1:
    #    event.time = dt.datetime.now()
    if postacq:
        post_acq = get_earliest_expected_dt(frame, event.time)
        if not post_acq:
            print('post event acquisition could not be identified')
            return False
        try:
            rc = lq.insert_new_eq2frame(eqid, fid, post_acq, active)
        except:
            print('error inserting to database')
    else:
        try:
            rc = lq.insert_new_eq2frame(eqid, fid, active = active)
        except:
            print('error inserting to database')
    if not rc:
        print('some problem importing event frame to eq2frame')
        return False
    return True


def is_blacklisted(eventid, blacklistfile = '/gws/nopw/j04/nceo_geohazards_vol1/public/LiCSAR_products/EQ/eq_blacklist.txt'):
    if os.path.exists(blacklistfile):
        if misc.grep1line(eventid,blacklistfile):
            print('the event {} is blacklisted'.format(eventid))
            return True
        else:
            return False
    else:
        print('WARNING: no blacklist file exists')
        return False


def process_all_eqs(minmag = 5.5, pastdays = 400, step = 2, overwrite = False):
    eventlist = get_eq_events(minmag, pastdays)
    print('identified {} events to process'.format(len(eventlist)))
    for event in eventlist:
        if not is_blacklisted(event.id):
            print(event.id)
            #keep them active if these are some late eqs..
            if event.time > dt.datetime.now()-timedelta(days=10):
                makeactive = True
            else:
                makeactive = False
            try:
                process_eq(event.id, step, overwrite, makeactive)
            except:
                print('some issue with event '+event.id)


def process_eq(eventid = 'us70008hvb', step = 1, overwrite = False, makeactive = False, skipchecks = False):
    event =  get_event_by_id(eventid)
    radius = get_range_from_magnitude(event.magnitude, event.depth, 'rad')
    if not radius:
        if skipchecks:
            print('the event is below threshold. setting minimal radius')
            radius = 0.18
        else:
            print('the event is below threshold. cancelling')
            return
    if not skipchecks:
        #check if the event is not blacklisted..
        blacklistfile = '/gws/nopw/j04/nceo_geohazards_vol1/public/LiCSAR_products/EQ/eq_blacklist.txt'
        if misc.grep1line(eventid,blacklistfile):
            print('the event is blacklisted. cancelling')
            return
    eventfile = os.path.join(public_path,'EQ',event.id+'.html')
    if not os.path.exists(eventfile):
        f=open(eventfile, "w")
        newline = "<a href='{0}'>USGS info on {1}</a> <br /> \n".format(event.url, event.id)
        f.write(newline)
        try:
            newline = "<a href='{0}'>USGS kml file</a> <br /> \n".format(event.getDetailURL().replace('geojson','kml'))
            f.write(newline)
        except:
            print('no KML exists for this event')
        f.close()
    frames = get_frames_in_event(event, radius)
    print('selected frames are:')
    print(frames)
    if len(frames) == 0:
        print('No frames are available for the event {0}'.format(event.id))
        return False
    eqid = import_to_licsinfo_eq(event, active = makeactive)
    for frame in frames:
        rc = import_to_licsinfo_eq2frame(eqid, event, frame[0], active = makeactive)
    if step == 1:
        print('{0} frames detected for event {1}, will be processing them'.format(str(len(frames)),event.id))
        indate = event.time-timedelta(days=max_days)
        offdate = event.time+timedelta(days=25)
        for frame in frames:
            frame = frame[0]
            track = str(int(frame[0:3]))
            if not os.path.exists(os.path.join(public_path,track,frame)):
                print('Frame '+frame+' was (probably) not initiated, trying to do it automatically')
                os.system('licsar_initiate_new_frame.sh {0} >/dev/null 2>/dev/null'.format(frame))
            framepubdir = os.path.join(public_path,track,frame)
            if os.path.exists(framepubdir):
                #if eqid:
                    #if equid means if first time to fille. then we should fill the eq2frame table
                #    try:
                #        import_to_licsinfo_eq2frame(eqid, event, frame)
                #    except:
                #        print('error during importing to eq2frame (database)')
                if not overwrite:
                    coseismic_ifgs = list_coseismic_ifgs(frame, event.time)
                    if not coseismic_ifgs:
                        start_processing=True
                    else:
                        print('there are some coseismic ifgs existing already, so skipping..')
                        start_processing=False
                else:
                    start_processing=True
                if start_processing:
                    print('..preparing frame {0} and sending processing jobs to LOTUS'.format(frame))
                    #print('licsar_make_frame.sh -n {0} 0 1 {1} {2}'.format('065D_05281_111313',str(indate.date()),str(offdate.date())))
                    #os.system('licsar_make_frame.sh -E -S {0} 1 1 {1} {2} >/dev/null 2>/dev/null'.format(frame,str(indate.date()),str(offdate.date())))
                    os.system('licsar_make_frame.sh -f -P -S {0} 1 1 {1} {2}'.format(frame,str(indate.date()),str(offdate.date())))
    elif step == 2:
        #second run - generate kmls:
        for frame in frames:
            frame = frame[0]
            track = str(int(frame[0:3]))
            #if not os.path.exists(os.path.join(procdir,'EQR',frame,'')
            new_kmls = create_kmls(frame,event.time, True, overwrite)
            if new_kmls:
                eqsfile = '/gws/nopw/j04/nceo_geohazards_vol1/public/LiCSAR_products/EQ/eqs.csv'
                if not misc.grep1line(event.id,eqsfile):
                    update_eq_csv(event.id, eqsfile)
                    update_eq2frames_csv(event.id)
                #update the event ID html file:
                print("Generated "+str(len(new_kmls))+" new (co-seismic ifg?) kml files")
                f=open(eventfile, "a+")
                for kml in new_kmls:
                    track = str(int(frame[0:3]))
                    fullwebpath = os.path.join(web_path, track, frame, 'interferograms', kml, frame+'_'+kml + '.kmz')
                    newline = '<a href="{0}">{1}: {2}.kmz</a> <br /> \n'.format(fullwebpath, frame, kml)
                    f.write(newline)
                    fullwebpath_metadata = os.path.join(web_path, track, frame, 'metadata')
                    newline = '<a href="{0}">{1} metadata</a> <br /> \n'.format(fullwebpath_metadata, frame)
                    f.write(newline)
                    fullwebpath_ifg = os.path.join(web_path, track, frame, 'interferograms', kml)
                    newline = '<a href="{0}">{1}: {2}</a> <br /> \n'.format(fullwebpath_ifg, frame, kml)
                    f.write(newline)
                f.close()
        if os.path.exists(eventfile):
            #this is to remove duplicities.. i know, should be done better..
            os.system('sort -u {0} > {0}.tmp; mv {0}.tmp {0}'.format(eventfile))
            print('done. Check this webpage:')
            print(os.path.join(web_path,'EQ',event.id+'.html'))


def main():
    #will process daily and check also for older earthquakes (to get max 12 days co-seismic pair)
    eventlist = get_eq_events(minmag) #search(starttime=datetime.now()-timedelta(days=max_days),
                           #endtime=datetime.now(),
                           #minmagnitude=5.5)
    #for debug only now
    #eventlist = search(starttime=datetime.now()-timedelta(days=5.2), endtime=datetime.now()-timedelta(days=4),minmagnitude=5.5)
    if eventlist:
        print('There were '+str(len(eventlist))+' events within last '+str(max_days)+' days')
        for event in eventlist:
            print('There was an event - ID '+event.id)
            #print(event.id)
            print("Time of event: "+str(event.time))
            print("Lat: "+str(event.latitude)+", Lon: "+str(event.longitude))
            print("Magnitude: "+str(event.magnitude)+", depth: "+str(event.depth)+" km")
            #print(event.location)
            #print(event.url)
            #now my functions:
            radius = get_range_from_magnitude(event.magnitude, event.depth, 'rad')
            #exceptions:
            if event.id in exceptions:
                print('ADDING EXCEPTION TO THIS EVENT')
                radius = get_range_from_magnitude(5.5, 10, 'rad')
            #print(radius)
            if radius:
                if event.hasProduct('shakemap'):
                    print('there is shakemap existing. we should download it and include to kml..')
                    print('see: '+event.url)
                    #print('detailURL is:'+event.getDetailURL())
                    #event.url is
                    #https://earthquake.usgs.gov/earthquakes/eventpage/us60004yps
                    #kml with eq info is:
                    #https://earthquake.usgs.gov/earthquakes/feed/v1.0/detail/us60004yps.kml
                    #event.getDetailURL() #here just change geojson to kml...
                    #event.url # shakemap can be downloaded with only contours > M5
                #global public_path
                #global web_path
                eventfile = os.path.join(public_path,'EQ',event.id+'.html')
                if not os.path.exists(eventfile):
                    f=open(eventfile, "w")
                    newline = "<a href='{0}'>USGS info on {1}</a> <br /> \n".format(event.url, event.id)
                    f.write(newline)
                    newline = "<a href='{0}'>USGS kml file</a> <br /> \n".format(event.getDetailURL().replace('geojson','kml'))
                    f.write(newline)
                    f.close()
                frames = get_frames_in_event(event, radius)
                #this would print info on next expected images:
                #get_next_expected_images(frames)
                if len(frames) == 0:
                    print('No frames are available for the event {0}'.format(event.id))
                else:
                    #ok, we will be processing this event, so let's:
                    #ingest it to the eq database
                    eqid = import_to_licsinfo_eq(event)
                    print('{0} frames detected for event {1}, will be processing them'.format(str(len(frames)),event.id))
                    for frame in frames:
                        frame = frame[0]
                        track = str(int(frame[0:3]))
                        indate = event.time-timedelta(days=max_days)
                        offdate = datetime.now()+timedelta(days=1)
                        if not os.path.exists(os.path.join(public_path,track,frame)):
                            print('Frame '+frame+' was (probably) not initiated, trying to do it automatically')
                            os.system('licsar_initiate_new_frame.sh {0} >/dev/null 2>/dev/null'.format(frame))
                        if os.path.exists(os.path.join(public_path,track,frame)):
                            if eqid:
                                #if equid means if first time to fille. then we should fill the eq2frame table
                                try:
                                    import_to_licsinfo_eq2frame(eqid, event, frame)
                                except:
                                    print('error during importing to eq2frame (database)')
                            print('..preparing frame {0} and sending processing jobs to LOTUS'.format(frame))
                            #print('licsar_make_frame.sh -n {0} 0 1 {1} {2}'.format('065D_05281_111313',str(indate.date()),str(offdate.date())))
                            os.system('licsar_make_frame.sh -E -P -S -N {0} 0 1 {1} {2} >/dev/null 2>/dev/null'.format(frame,str(indate.date()),str(offdate.date())))
                            #this way is a kind of prototype and should be improved
                            new_kmls = create_kmls(frame,event.time)
                            #try:
                            #    new_kmls = create_kmls(frame,event.time)
                            #except:
                            #    new_kmls = None
                            #    print('some error happened during processing frame '+frame)
                            if new_kmls:
                                #so if something was generated, update it to the list of eqs
                                eqsfile = '/gws/nopw/j04/nceo_geohazards_vol1/public/LiCSAR_products/EQ/eqs.csv'
                                if not misc.grep1line(event.id,eqsfile):
                                    update_eq_csv(event.id, eqsfile)
                                    update_eq2frames_csv(event.id)
                                #update the event ID file:
                                print("Generated "+str(len(new_kmls))+" new (co-seismic ifg?) kml files")
                                f=open(eventfile, "a+")
                                for kml in new_kmls:
                                    track = str(int(frame[0:3]))
                                    fullwebpath = os.path.join(web_path, track, frame, 'interferograms', kml, frame+'_'+kml + '.kmz')
                                    newline = '<a href="{0}">{1}: {2}.kmz</a> <br /> \n'.format(fullwebpath, frame, kml)
                                    f.write(newline)
                                    fullwebpath_metadata = os.path.join(web_path, track, frame, 'metadata')
                                    newline = '<a href="{0}">{1} metadata</a> <br /> \n'.format(fullwebpath_metadata, frame)
                                    f.write(newline)
                                    fullwebpath_ifg = os.path.join(web_path, track, frame, 'interferograms', kml)
                                    newline = '<a href="{0}">{1}: {2}</a> <br /> \n'.format(fullwebpath_ifg, frame, kml)
                                    f.write(newline)
                                f.close()
                        else:
                            print('The frame '+frame+' has not been initialized properly using automatic approach')
                    if os.path.exists(eventfile):
                        os.system('sort -u {0} > ~/tmpbardak; mv ~/tmpbardak {0}'.format(eventfile))
                        print('done. Check this webpage:')
                        print(os.path.join(web_path,'EQ',event.id+'.html'))

if __name__ == '__main__':
    main()
