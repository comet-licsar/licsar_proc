#!/usr/bin/env python
import pandas as pd
import sys, os
import datetime as dt
#import framecare as fc

currentdir = os.environ['LiCSAR_procdir']
pubdir = os.environ['LiCSAR_public']

frame=sys.argv[1]
try:
    infile=sys.argv[2]
except:
    #let's try getting gaps
    track = str(int(frame[0:3]))
    pubpath = os.path.join(pubdir,track,frame)
    networkfile = os.path.join(pubpath,'metadata','network.png')
    gapfile = os.path.join(pubpath,'metadata','gaps.txt')
    print('updating network plot and gap list')
    os.system('module load LiCSBAS; plot_network.py {0} {1} {2}'.format(pubpath,networkfile,gapfile))
    if not os.path.exists(gapfile):
        print('error generating gapfile, exiting')
        exit()
    else:
        print('gapfile generated')
        infile = gapfile
        splitter = '_'

with open(infile, 'r') as f:
    x = f.readlines()
    #from framecare import get_ifg_list_pubdir
    #x = get_ifg_list_pubdir(framename)
    #splitter = '_'

epoch1 = []
epoch2 = []
for a in x:
    if '20' in a:
        epoch1.append(int(a.split(splitter)[0]))
        epoch2.append(int(a.split(splitter)[1]))
    if 'Connected network ' in a:
        break
#tolerance of 25 days should include 2-4 epochs
tol = pd.Timedelta('25 days')
#licsar_params = '-S'
#licsar_toggles = '1 1'

gaps = pd.DataFrame()
gaps['epoch1'] = pd.to_datetime(epoch1,format='%Y%m%d') - tol
gaps['epoch2'] = pd.to_datetime(epoch2,format='%Y%m%d') + tol
gaps_corrected = gaps.copy()
gapscorindex = 0
#print(gaps)
for x in range(len(gaps)-1):
    if gaps.epoch2[x] >= gaps.epoch1[x+1]:
        gaps_corrected.epoch2[x - gapscorindex] = gaps.epoch2[x+1]
        gaps_corrected = gaps_corrected.drop(x+1)
        gapscorindex = gapscorindex + 1


gaps_corrected = gaps_corrected.reset_index(drop=True)
gaps_corrected['btemp'] = gaps_corrected.epoch2 - gaps_corrected.epoch1
print('gap dates in the frame')
print('(frame,datein,dateout,total btemp)')
for x in range(len(gaps_corrected)):
    date1 = str(gaps_corrected.epoch1[x].date())
    date2 = str(gaps_corrected.epoch2[x].date())
    btemp = str(gaps_corrected.btemp[x].days)
    #print('licsar_make_frame.sh {0} {1} {2} {3} {4}'.format(
    #    licsar_params,frame,licsar_toggles,date1,date2))
    print('{0},{1},{2},{3}'.format(
        frame,date1,date2,btemp))
