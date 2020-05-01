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

public_path = os.environ['LiCSAR_public']
procdir_path = os.environ['LiCSAR_procdir']
web_path = 'http://gws-access.ceda.ac.uk/public/nceo_geohazards/LiCSAR_products/'

#you may want to change these parameters:
max_days = 30
minmag = 5.5
#here will be exceptions (i.e. eqs that MUST be processed):
#minmag = 4.6
exceptions = ['us60008e8e', 'us70008dx7']


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
    eq_df['link'] = "<a href='http://gws-access.ceda.ac.uk/public/nceo_geohazards/LiCSAR_products/EQ/" + eq_df['USGS_ID'] + ".html' target='_blank'>Link<a>"
    dbcols = ['USGS_ID','magnitude','depth','time','lat','lon', 'link']
    eq_df[dbcols].to_csv(csvfile, sep = ';', index=False)
    return True

def update_eq_csv(eventid, csvfile = '/gws/nopw/j04/nceo_geohazards_vol1/public/LiCSAR_products/EQ/eqs.csv'):
    #get all eqs
    query = "select * from eq where USGS_ID='{}';".format(eventid)
    eq_df = lq.do_pd_query(query)
    eq_df['link'] = "<a href='http://gws-access.ceda.ac.uk/public/nceo_geohazards/LiCSAR_products/EQ/" + eq_df['USGS_ID'] + ".html' target='_blank'>Link<a>"
    dbcols = ['USGS_ID','magnitude','depth','time','lat','lon', 'link']
    eq_df[dbcols].to_csv(csvfile, mode='a', sep = ';', index = False, header=False)
    return True

def update_eq2frames_csv(eventid, csvfile = '/gws/nopw/j04/nceo_geohazards_vol1/public/LiCSAR_products/EQ/eqframes.csv'):
    """
       This csv will be loaded to the eq responder map
    """
    #fid = lq.get_frame_polyid(framename)
    #try:
    #    fid = lq.sqlout2list(fid)[0]
    #except:
    #    print('error - the frame {} does not exist in licsinfo database'.format(framename))
    #query = "select pg.geom from eq2frame e2f inner join eq on e2f.eqid=eq.eqid inner join polygs2gis pg on pg.polyid=e2f.fid where eq.USGS_ID='{0}' and e2f.fid={1};".format(eventid, str(fid))
    #query = "select pg.geom,e2f.fid from eq2frame e2f inner join eq on e2f.eqid=eq.eqid inner join polygs2gis pg on pg.polyid=e2f.fid where eq.USGS_ID='{0}'".format(eventid)
    query = "select p.polyid_name as frame, aswkb(pg.geom) as the_geom, eq.USGS_ID as usgsid, eq.location from eq2frame e2f inner join eq on e2f.eqid=eq.eqid inner join polygs2gis pg on pg.polyid=e2f.fid inner join polygs p on p.polyid=e2f.fid where eq.USGS_ID='{0}'".format(eventid)
    #hmm.. this does not load geometry as wkb... need to solve it!
    e2f = lq.do_pd_query(query)
    if len(e2f) < 1:
        print('error - nothing found in eq2frames for this event {}'.format(eventid))
        return False
    e2f['track'] = e2f.frame.apply(lambda x:str(int(x[:3])))
    e2f['download'] = e2f['track']*0
    for i, f in e2f.iterrows():
        try:
            row = fc.frame2geopandas(f['frame']).iloc[0]
            geometry = row['geometry'].wkb_hex
            e2f.iloc[i].the_geom = geometry
        except:
            print('probably this frame is not in geom database')
            return False
        e2f.iloc[i]['download'] = "<a href='http://gws-access.ceda.ac.uk/public/nceo_geohazards/LiCSAR_products/{0}/{1} target='_blank'>Link<a>".format(f['track'], f['frame'])
    #e2f['download'] = "<a href='http://gws-access.ceda.ac.uk/public/nceo_geohazards/LiCSAR_products/{0}/{1} target='_blank'>Link<a>".format(track, framename)
    dbcols = ['the_geom','frame','usgsid','download', 'location']
    e2f[dbcols].to_csv(csvfile, mode='a', sep = ';', index = False, header=False)
    return True

def create_kmls(frame, toi):
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
            filelist = lq.get_frame_files_date(frame, date)
            tof = lq.get_time_of_file(filelist[0][1])
            if tof > toi:
                is_coseismic = True
            if tof < toi:
                is_preseismic = True
        if is_coseismic or is_postseismic:
            if not os.path.exists(os.path.join(products_path,pifg,frame+'_'+pifg+'.kmz')):
                #first of all regenerate preview in full resolution - and include unfiltered data
                os.system('create_geoctiffs_to_pub.sh -F -u {0} {1}'.format(frameproc_path, pifg))
                print('to generate hires geotiffs:')
                print('create_geoctiffs_to_pub.sh -H -u {0} {1}'.format(frameproc_path, pifg))
                print('cd {0}/GEOC_50m/{1}; for x in \`ls *tif\`; do mv \$x \$LiCSAR_public/.../{1}/\`echo \$x | sed \'s/tif/hires.tif/\'\`; cd -'.format(frameproc_path, pifg))
                #print('mv {0}/{1}/*tif')
                #generate kmz
                #os.system('create_kmz.sh {0}'.format(os.path.join(products_path,pifg))) #this will do it directly in LiCSAR_public
                os.system('create_kmz.sh {0}'.format(os.path.join(frameproc_path,'GEOC',pifg)))
                #if os.path.exists(os.path.join(products_path,pifg,pifg+'.kmz')):
                if os.path.exists(os.path.join(frameproc_path,'GEOC',pifg,pifg+'.kmz')):
                    #os.rename(os.path.join(products_path,pifg,pifg+'.kmz'), os.path.join(products_path,pifg,frame+'_'+pifg+'.kmz'))
                    os.rename(os.path.join(frameproc_path,'GEOC',pifg,pifg+'.kmz'), os.path.join(products_path,pifg,frame+'_'+pifg+'.kmz'))
                    #clean the generated hires files
                    shutil.rmtree(os.path.join(frameproc_path,'GEOC',pifg))
                    selected_ifgs.append(pifg)
    return selected_ifgs

def get_earliest_expected_dt(frame, eventtime):
    masterdate = fc.get_master(frame, asdatetime = True)
    if not masterdate:
        print('error getting masterdate')
        return False
    daysdiff = (eventtime - masterdate).days
    noepochs = int(np.floor(daysdiff/6))
    lastepoch = masterdate+dt.timedelta(days=noepochs*6)
    expected_dt = lastepoch+dt.timedelta(days=6)
    return expected_dt

def get_next_expected_datetime(frame, eventtime):
    track = str(int(frame[0:3]))
    lastmonth = eventtime-timedelta(days=30)
    eqimages = get_images_for_frame(frame, startdate = lastmonth.date(), enddate = eventtime.date(), sensType = 'IW')
    eqimgdates = set()
    for eqi in eqimages:
        imgdate = eqi.split('T')[0].split('_')[-1]
        imgdate = datetime.strptime(imgdate,'%Y%m%d')
        eqimgdates.add(imgdate)
    eqimgdates = sorted(eqimgdates) # list now
    
    #lastimage = eqimgdates[-1]
    expectedtime = eqimages[-1].split('T')[1].split('_')[0]
    expectedtime = datetime.strptime(expectedtime[0:4],'%H%M')
    nextposimage = eqimgdates[-1]+timedelta(days=6)
    nextposimage = nextposimage.replace(hour=expectedtime.hour, minute=expectedtime.minute)
    
    imgdiffs = []
    for i in range(len(eqimgdates)-1):
        imgdiffs.append((eqimgdates[i+1] - eqimgdates[i]).days)
    bestcasediff = min(imgdiffs)
    worstcasediff = max(imgdiffs)
    if not worstcasediff == bestcasediff:
        print('Temporal sampling is irregular (expecting worst case), i.e. '+str(worstcasediff)+' days diff')
    imgdiff = worstcasediff
    
    nextexpectedimage = eqimgdates[-1] + timedelta(days=imgdiff)
    nextexpectedimage = nextexpectedimage.replace(hour=expectedtime.hour, minute=expectedtime.minute)
    
    returnlist = [nextexpectedimage, nextposimage]
    return returnlist

def get_next_expected_images(frames, eventtime):
    print('checking first expected post-seismic acquisition time for covering frames:')
    for frame in frames:
        frame = frame[0]
        print(frame+': ')
        [nextexp, nextpos] = get_next_expected_datetime(frame, eventtime)
        print('Expected: '+str(nextexp))
        if nextexp != nextpos:
            print('Warning, this area is not observed at the highest frequency')
            print('Next possible flight: '+str(nextpos))

def get_eq_events(minmag = 5.5, max_days = 30):
    return search(starttime=datetime.now()-timedelta(days=max_days),
                           endtime=datetime.now(),
                           minmagnitude=minmag)

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

def import_to_licsinfo_eq(event):
    if lq.get_eqid(event.id):
        print('eq already in database, cancelling')
        return False
    else:
        eqid = lq.insert_new_eq(event)
        return eqid

def import_to_licsinfo_eq2frame(eqid, event, frame, postacq = True):
    fid = lq.get_frame_polyid(frame)
    try:
        fid = lq.sqlout2list(fid)[0]
    except:
        print('the frame does not exist in licsinfo!')
        return False
    if postacq:
        post_acq = get_earliest_expected_dt(frame, event.time)
        if not post_acq:
            print('post event acquisition could not be identified')
            return False
        rc = lq.insert_new_eq2frame(eqid, fid, post_acq)
    else:
        rc = lq.insert_new_eq2frame(eqid, fid)
    if not rc:
        print('some problem importing event frame to eq2frame')
        return False
    return True

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
                            os.system('licsar_make_frame.sh -S -N {0} 0 1 {1} {2} >/dev/null 2>/dev/null'.format(frame,str(indate.date()),str(offdate.date())))
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
