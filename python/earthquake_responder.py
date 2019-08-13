#!/usr/bin/env python

import os
import LiCSquery as lq
import datetime
from datetime import datetime, timedelta
from libcomcat.search import search,count,get_event_by_id
import pandas as pd

public_path = os.environ['LiCSAR_public']
web_path = 'http://gws-access.ceda.ac.uk/public/nceo_geohazards/LiCSAR_products/'

def get_range_from_magnitude(M, depth, unit = 'km'):
    #M = 5.5
    #depth = 9 #km
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
    if M > 8.5 and depth <= 250:
        distance = 2500
    else:
        distance = eq_limits.query('magnitude == {0} and depth >= {1}'.format(M,depth))['distance']
        if len(distance)>0:
            distance = int(distance.to_string().split()[1])
        else:
            distance = None
            print('The earthquake parameters do not fit with the limit table')
    #rad or km?
    if unit == 'rad':
        distance = (360.0 / 40007.86) * distance
    return distance

def create_kmls(frame, toi):
    #toi is Time Of Interest - and should be as datetime
    doi = toi.date()
    track = str(int(frame[0:3]))
    global public_path
    products_path = os.path.join(public_path, track, frame, 'products')
    doi_str = int(doi.strftime('%Y%m%d'))
    publicifgs = os.listdir(products_path)
    ## should do check for valid folder name here!
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
            tof = get_time_of_file(filelist[0][1])
            if tof < toi:
                is_coseismic = True
        if (doi_str == slv):
            date = datetime.strptime(str(slv),'%Y%m%d').date()
            filelist = lq.get_frame_files_date(frame, date)
            tof = lq.get_time_of_file(filelist[0][1])
            if tof > toi:
                is_coseismic = True
        if is_coseismic:
            if not os.path.exists(os.path.join(products_path,pifg,pifg+'.kmz')):
                #generate kmz
                os.system('create_kmz.sh {0}'.format(os.path.join(products_path,pifg)))
                selected_ifgs.append(pifg)
    return selected_ifgs

def main():
    #will process daily and check also for older earthquakes (to get max 12 days co-seismic pair)
    eventlist = search(starttime=datetime.now()-timedelta(days=13),
                           endtime=datetime.now(),
                           minmagnitude=5.5)
    #for debug only now
    eventlist = search(starttime=datetime.now()-timedelta(days=5.2), endtime=datetime.now()-timedelta(days=4),minmagnitude=5.5)
    if eventlist:
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
            #print(radius)
            if event.hasProduct('shakemap'):
                print('there is shakemap existing. we should download it and include to kml..')
                print('see: '+event.url)
                #event.getDetailURL() #here just change geojson to kml...
                #event.url # shakemap can be downloaded with only contours > M5
            minlon = event.longitude - radius
            maxlat = event.latitude + radius
            maxlon = event.longitude + radius
            minlat = event.latitude - radius
            #global public_path
            #global web_path
            eventfile = os.path.join(public_path,'EQ',event.id+'.html')
            frames = lq.get_frames_in_polygon(minlon,maxlon,minlat,maxlat)
            print('{0} frames detected for event {1}, will be processing them'.format(str(len(frames)),event.id))
            for frame in frames:
                frame = frame[0]
                indate = event.time-timedelta(days=13)
                offdate = datetime.now()
                print('..preparing frame {0} and sending processing jobs to LOTUS'.format(frame))
                #print('licsar_make_frame.sh -n {0} 0 1 {1} {2}'.format('065D_05281_111313',str(indate.date()),str(offdate.date())))
                os.system('licsar_make_frame.sh -S -N {0} 0 1 {1} {2} >/dev/null 2>/dev/null'.format(frame,str(indate.date()),str(offdate.date())))
                #this way is a kind of prototype and should be improved
                new_kmls = create_kmls(frame,event.time)
                if new_kmls:
                    #update the event ID file:
                    print("Generated "+str(len(new_kmls))+" new (co-seismic ifg?) kml files")
                    f=open(eventfile, "a+")
                    for kml in new_kmls:
                        track = str(int(frame[0:3]))
                        fullwebpath = os.path.join(web_path, track, frame, 'products', kml, kml + '.kmz')
                        newline = '<a href="{0}">{1}: {2}</a> <br /> \n'.format(fullwebpath, frame, kml)
                        f.write(newline)
                    f.close()
            if os.path.exists(eventfile):
                os.system('sort -u {0} > ~/tmpbardak; mv ~/tmpbardak {0}'.format(eventfile))
                print('done. Check this webpage:')
                print(os.path.join(web_path,'EQ',event.id+'.html'))

if __name__ == '__main__':
    main()
