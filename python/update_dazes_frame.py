#!/usr/bin/env python3

# import framecare as fc
import daz_lib_licsar as dll
import os
import pandas as pd


if len(sys.argv) < 2:
    print('Please provide frame id: i.e python update_dazes_frame.py 021D_05266_252525')
    sys.exit(1)

##variables
frame=sys.argv[1]

tr = int(frame[:3])
metadir = os.path.join(os.environ['LiCSAR_public'], str(tr), frame, 'metadata')
dazcsv = os.path.join(metadir, frame+'.dazdrg')
if not os.path.exists(metadir):
    print('metadata folder does not exist - cancelling')
    exit(1)

# get all
dv=dll.get_daz_frame(frame) #, include_corrections=True)
dv=dv.drop('polyid', axis=1)
dv['epoch'] = dv.apply(lambda x: pd.to_datetime(x['epoch']).date(), axis=1)

if os.path.exists(dazcsv):
    dvcsv = pd.read_csv(dazcsv)
    dvcsv['epoch'] = dvcsv.apply(lambda x: pd.to_datetime(x['epoch']).date(), axis=1)
    for exepoch in dvcsv.epoch.values:
        dv = dv[dv.epoch != exepoch]

if dv.empty:
    print('dazdrg file is up to date - nothing to change')
    exit()

dvout = pd.DataFrame()
for epoch in dv.epoch.values:
    # update it but only if we have gacos data:
    epochstr = epoch.strftime('%Y%m%d')
    gacosfile = os.path.join(os.environ['LiCSAR_public'], str(tr), frame, 'epochs', epochstr, epochstr+'.sltd.geo.tif')
    if not os.path.exists(gacosfile):
        print('no gacos for '+epochstr+' - skipping')
        continue
    dvepoch = dll.get_daz_frame(
            frame,
            fulloutput=True,
            include_corrections=True,
            use_iri_hei=False,
            corr_per_swath=True,
            datemin=epoch,
            datemax=epoch
        )
    dvout = pd.concat([dvout, dvepoch])

dv = dvout.reset_index(drop=True)
dv['daz_mm'] = dv['daz']*14000
dv['daz_icc_mm'] = dv['cc_azi']*14000
dv['drg_icc_mm'] = dv['cc_range']*2300
dv['daz_iono_mm'] = dv['daz_iono']*14000
dv['daz_SET_mm'] = dv['daz_SET']*14000
dv = dv.drop(['daz', 'cc_azi', 'cc_range', 'daz_iono', 'daz_SET'], axis=1)
dv=dv.drop('polyid', axis=1)
dv=pd.concat([dvcsv, dv])
dv.to_csv(dazcsv, index=False, float_format='%.2f')




