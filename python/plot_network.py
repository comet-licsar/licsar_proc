#!/usr/bin/env python3
# two outputs are: png plot and text file for gaps (both are mandatory)
# e.g. $LiCSAR_public/$track/$frame output.png gaps.txt
# then... there is optional parameter for check_common_bursts --- so if you add '1' , it will additionally check for the common bursts. will take more time though...
#%% Import
import os
import sys
import numpy as np
import LiCSBAS_io_lib as io_lib
#import LiCSBAS_plot_lib as plot_lib
from LiCSBAS_plot_lib import *
import LiCSBAS_tools_lib as tools_lib
import LiCSBAS_inv_lib as inv_lib
import datetime as dt
import s1data as s1
import pandas as pd
import framecare as fc

check_common_bursts = False

#%%
def read_bperp_file(bperp_file, imdates, return_missflag = False):
    """
    updated from LiCSBAS io_lib function to give 0 for missing imdates
    
    bperp_file (baselines) contains (m: primary (master), s: secondary,
                                     sm: single prime):
          smdate    sdate    bp    dt
        20170302 20170326 130.9  24.0
        20170302 20170314  32.4  12.0
    Old bperp_file contains (m: primary (master), s:secondary,
                             sm: single prime):
        num    mdate    sdate   bp   dt  dt_m_sm dt_s_sm bp_m_sm bp_s_sm
          1 20170218 20170326 96.6 36.0    -12.0    24.0    34.2   130.9
          2 20170302 20170314 32.4 12.0      0.0    12.0     0.0    32.4
    Return: bperp
    """
    bperp = []
    missflag = False
    bperp_dict = {}
    ### Determine type of bperp_file; old or not
    with open(bperp_file) as f:
        line = f.readline().split() #list
    if len(line) == 4: ## new format
        bperp_dict[line[0]] = '0.00' ## single prime. unnecessary?
        with open(bperp_file) as f:
            for l in f:
                if len(l.split()) == 4:
                    bperp_dict[l.split()[1]] = l.split()[2]
    else: ## old format
        with open(bperp_file) as f:
            for l in f:
                bperp_dict[l.split()[1]] = l.split()[-2]
                bperp_dict[l.split()[2]] = l.split()[-1]
    for imd in imdates:
        if imd in bperp_dict:
            bperp.append(float(bperp_dict[imd]))
        else: ## If no key exists
            bperp.append(0)
            missflag = True
            if not return_missflag:
                print('WARNING: bperp for {} not found, nullifying'.format(imd))
            #return False
    if return_missflag:
        return bperp, missflag
    else:
        return bperp



#%%
def plot_network_upd(ifgdates, bperp, frame, pngfile, firstdate = dt.datetime(2014, 9, 25), lastdate = dt.datetime(2025, 12, 31)):
    """
    Plot network of interferometric pairs.
    bperp can be dummy (-1~1).
    Suffix of pngfile can be png, ps, pdf, or svg.
    """
    imdates_all = tools_lib.ifgdates2imdates(ifgdates)
    n_im_all = len(imdates_all)
    imdates_dt_all = np.array(([dt.datetime.strptime(imd, '%Y%m%d') for imd in imdates_all])) ##datetime
    ifgdates = list(set(ifgdates))
    ifgdates.sort()
    imdates = tools_lib.ifgdates2imdates(ifgdates)
    n_im = len(imdates)
    imdates_dt = np.array(([dt.datetime.strptime(imd, '%Y%m%d') for imd in imdates])) ##datetime
    #
    ### Identify gaps    
    G = inv_lib.make_sb_matrix(ifgdates)
    ixs_inc_gap = np.where(G.sum(axis=0)==0)[0]
    #
    ### Plot fig
    #figsize_x = np.round(((imdates_dt_all[-1]-imdates_dt_all[0]).days)/80)+2
    figsize_x = np.round(((lastdate-firstdate).days)/80)+2
    #fig = plt.figure(figsize=(figsize_x, 6))
    fig = plt.figure(figsize=(figsize_x, 7))
    #ax = fig.add_axes([0.06, 0.12, 0.92,0.85])
    ax = fig.add_axes([0.03, 0.12, 0.94,0.8])
    #
    ### IFG blue lines
    for i, ifgd in enumerate(ifgdates):
        ix_m = imdates_all.index(ifgd[:8])
        ix_s = imdates_all.index(ifgd[-8:])
        label = 'IFG' if i==0 else '' #label only first
        plt.plot([imdates_dt_all[ix_m], imdates_dt_all[ix_s]], [bperp[ix_m],
                bperp[ix_s]], color='b', alpha=0.6, zorder=2, label=label)
    #
    #
    ### Image points and dates
    ax.scatter(imdates_dt_all, bperp, alpha=0.6, zorder=4)
    for i in range(n_im_all):
        if bperp[i] > np.median(bperp): va='bottom'
        else: va = 'top'
        ax.annotate(imdates_all[i][4:6]+'/'+imdates_all[i][6:],
                    (imdates_dt_all[i], bperp[i]), ha='center', va=va, zorder=8)
    #
    #
    ### gaps
    if len(ixs_inc_gap)!=0:
        gap_dates_dt = []
        for ix_gap in ixs_inc_gap:
            ddays_td = imdates_dt[ix_gap+1]-imdates_dt[ix_gap]
            gap_dates_dt.append(imdates_dt[ix_gap]+ddays_td/2)
        plt.vlines(gap_dates_dt, 0, 1, transform=ax.get_xaxis_transform(),
                   zorder=1, label='Gap', alpha=0.6, colors='k', linewidth=3)
    #
    #
    ### Locater        
    loc = ax.xaxis.set_major_locator(mdates.AutoDateLocator())
    try:  # Only support from Matplotlib 3.1
        ax.xaxis.set_major_formatter(mdates.ConciseDateFormatter(loc))
    except:
        ax.xaxis.set_major_formatter(mdates.DateFormatter('%Y/%m/%d'))
        for label in ax.get_xticklabels():
            label.set_rotation(20)
            label.set_horizontalalignment('right')
    #
    #
    ax.grid(b=True, which='major')
    ### Add bold line every 1yr
    ax.xaxis.set_minor_locator(mdates.YearLocator())
    ax.grid(b=True, which='minor', linewidth=2)
    ax.set_xlim((firstdate, lastdate))
    #ax.set_xlim((imdates_dt_all[0]-dt.timedelta(days=10),
    #             imdates_dt_all[-1]+dt.timedelta(days=10)))
    ### Labels and legend
    plt.xlabel('Time')
    if np.all(np.abs(np.array(bperp))<=1): ## dummy
        plt.ylabel('dummy')
    else:
        plt.ylabel('Bperp [m]')
    #
    # 2022-04-19 adding dots of 'existing epochs'
    print("getting existing epochs for the frame bounding box")
    epochdates = s1.get_epochs_for_frame(frame, firstdate.date(), lastdate.date(), returnAsDate = True)
    epochdates.sort()
    for imd in imdates_dt:
        imdd = imd.date()
        if imdd in epochdates:
            #print('debug - found and removed ok: '+str(imdd))
            epochdates.remove(imdd)
    if check_common_bursts:
        print("checking if the epochs have any common burst - if not, will plot them anyway, in gray")
        framebursts = fc.lq.sqlout2list(fc.get_bidtanxs_in_frame(frame))
        epochdates_outburst = []
        for imdd in epochdates:
            if len(fc.get_frame_files_date(frame, imdd))>0:
                continue  # if we find something in database, it just means there was this acquisition, so plot it
            print('checking epoch '+str(imdd))
            # there is some overlap but does it have the same bursts?
            try:
                images = s1.get_images_for_frame(frame, startdate = imdd-dt.timedelta(days=1), enddate = imdd+dt.timedelta(days=1), asf = False)
                for im in images:
                    bursts = fc.lq.sqlout2list(fc.get_bursts_in_file(im))
                    if not bursts:
                        filepath = s1.get_neodc_path_images(im, file_or_meta=True)[0]
                        _ = fc.ingest_file_to_licsinfo(filepath, isfullpath=True)
                        bursts = fc.lq.sqlout2list(fc.get_bursts_in_file(im))
                    isinframe = False
                    for b in bursts:
                        if b in framebursts:
                            isinframe = True
                            break
                    if not isinframe:
                        epochdates.remove(imdd)  # remove from getting plotted as red circles
                        epochdates_outburst.append(imdd)
            except:
                print('some error double checking epoch '+str(imdd)+'. keeping it')
        if epochdates_outburst:
            ax.scatter(epochdates_outburst, np.zeros(len(epochdates_outburst)), facecolors='none', edgecolors='gray', label='existing acquisition with no frame burst (debug)')
    if epochdates:
        ax.scatter(epochdates,np.zeros(len(epochdates)), facecolors='none', edgecolors='red', label='existing acquisition')
    # adding timestamp
    timestamp = 'updated: '+str(dt.datetime.now().strftime("%Y-%m-%d %I:%M:%S"))
    plt.title(frame+', '+timestamp)
    ax.title.set_size(16)
    #plt.text(0.5,0.5,timestamp)
    plt.legend()
    ### Save
    plt.savefig(pngfile)
    plt.close()





#%% File setting
try:
    framedir = sys.argv[1]
    pngfile = sys.argv[2] 
    gapfile = sys.argv[3]
    try:
        islast = sys.argv[4]
        if islast == '1':
            check_common_bursts = True
            print('will check for the common bursts (more correct "red circles")')
    except:
        pass
except:
    print('Usage: ')
    print('plot_network.py path_to_frame_directory out_png_file out_gaps_file [1]')
    print('(where the optional 1 will mean (long but correct) checking of the frame-related epochs based on bursts)')
    exit()


ifgdir = os.path.join(framedir, 'interferograms')
bperp_file = os.path.join(framedir, 'metadata', 'baselines')
if not os.path.exists(ifgdir):
    # update to have it work in BATCH_CACHE_DIR
    ifgdir = os.path.join(framedir, 'GEOC')
    bperp_file = os.path.join(framedir, 'baselines')
    if not os.path.exists(ifgdir):
        print('error, no interferograms found for this frame')
        exit()
    else:
        print('generating baselines file in custom directory')
        cmd = 'cd {0}; mk_bperp_file.sh; mv baselines {1} 2>/dev/null'.format(framedir, bperp_file)
        rc = os.system(cmd)


if os.path.exists(pngfile):
    os.remove(pngfile)
if os.path.exists(gapfile):
    os.remove(gapfile)
#%%
ifgdates = tools_lib.get_ifgdates(ifgdir)
imdates = tools_lib.ifgdates2imdates(ifgdates)


if not os.path.exists(bperp_file):
    print('No baselines file exists. The Bperps will be estimated')
    frame = os.path.basename(framedir)
    bpd = fc.make_bperp_file(frame, bperp_file, donotstore=False)

# else:
# horrible fix but seems necessary...
rc = os.system("sed -i 's/\.0//g' "+bperp_file)
#    print('Make dummy bperp')
#    bperp_file = os.path.join(framedir,'baselines_tmp.txt')
#    io_lib.make_dummy_bperp(bperp_file, imdates)



try:
    bperp, ismissing = read_bperp_file(bperp_file, imdates, return_missflag = True)
    try:
        # just in case...
        bperp = np.array(bperp)
        bperp[np.isnan(bperp)] = 0
        absbp = np.abs(bperp)
        over = np.where(absbp>800)
        if len(over)>0:
            print('WARNING, removing bperps that are over threshold of 800 m. These will be reestimated.')
            bperp[over]=0
            ismissing = True
    except:
        print('an error trying to remove bperps over 800 m')

    # double check missing - count zeroes
    if not ismissing:
        if len(bperp)>1:
            absbp=np.abs(bperp)
            absbp.sort()
            if absbp[1] == 0:
                ismissing = True
    if ismissing:
        print('some epochs have missing bperps, trying to find them through ASF')
        frame=os.path.basename(framedir)
        bperp=np.array(bperp)
        imdates=np.array(imdates)
        missingdates = imdates[bperp==0]
        missingdates2 = imdates[np.abs(bperp)>400]
        missingdates = np.concatenate((missingdates2, missingdates))
        missingdates2 = imdates[np.isnan(bperp)]
        missingdates = np.concatenate((missingdates2, missingdates))
        refdate = fc.get_master(frame)
        missingdates = missingdates[missingdates != refdate]
        # load existing
        prevbp = pd.read_csv(bperp_file, header=None, sep = ' ')
        prevbp.columns = ['ref_date', 'date', 'bperp', 'btemp']
        #print('TODO - remove missingepochs from prevbp')
        # get new - try first only from ASF (more accurate)
        bpd = fc.make_bperp_file(frame, bperp_file, asfonly = True, donotstore = True)
        stillmissing = []
        for m in missingdates:
            # first drop it from the prevbp:
            mint = int(m)
            prevbp = prevbp.drop(prevbp[prevbp.date == mint].index)
            mpd = bpd[bpd.date==m]
            if not mpd.empty:
                mbperp = mpd.bperp.mean()
                mbtemp = mpd.btemp.values[0]
                prevbp.loc[len(prevbp.index)] = [int(refdate), int(m), mbperp, int(mbtemp)] # new line
            else:
                #mbperp = 0
                #mbtemp = fc.datediff(refdate, m)
                print('no ASF information for epoch '+m+'. Adding for LiCSAR estimation.') #Storing only bperp=0')
                stillmissing.append(m)
                ''' NOT COMPLETE YET - SOMETHING IS WRONG IN THIS BELOW:
                print('no ASF information for epoch '+m+'. Estimating from LiCSAR db - slow way now') #Storing only bperp=0')
                try:
                    mepl, mbperpl = fc.get_bperp_estimates(frame, epochs = [m])
                    mbperp = round(mbperpl[0])
                except:
                    print('ERROR for epoch '+m+'. Setting zero.')
                    mbperp = 0
                '''
        if stillmissing:
            bperps = fc.estimate_bperps(frame, stillmissing, return_epochsdt=False)
            bperps = np.array(bperps).astype(int)
            i = 0
            for m in stillmissing:
                mbperp = bperps[i]
                mbtemp = fc.datediff(refdate, m)
                prevbp.loc[len(prevbp.index)] = [int(refdate), int(m), int(mbperp), int(mbtemp) ]
                i = i+1
        prevbp = prevbp.sort_values('btemp').reset_index(drop=True)
        #bpd.to_csv(bperp_file, sep = ' ', index = False, header = False)
        # bperps = bperps.astype(int)  # for some reason we still export as floats!
        for col in prevbp.columns:
            prevbp[col] = prevbp[col].astype(int)
        prevbp.to_csv(bperp_file, sep = ' ', index = False, header = False)
        bperp = read_bperp_file(bperp_file, imdates)
except:
    print('error reading baselines file! trying to fully recreate through ASF')
    try:
        if os.path.exists(bperp_file):
            os.remove(bperp_file)
        frame=os.path.basename(framedir)
        rc = fc.make_bperp_file(frame, bperp_file)
        bperp = read_bperp_file(bperp_file, imdates)
    except:
        print('some error occurred. Making dummy bperp')
        bperp_file = os.path.join(framedir,'baselines_tmp.txt')
        io_lib.make_dummy_bperp(bperp_file, imdates)
        bperp = read_bperp_file(bperp_file, imdates)


frame = os.path.basename(framedir)
plot_network_upd(ifgdates, bperp, frame, pngfile)
os.system('chmod 777 '+pngfile+' 2>/dev/null')
os.system('chmod 777 '+bperp_file+' 2>/dev/null')
rc = os.system("sed -i 's/\.0//g' "+bperp_file)  # just in case...

## Identify gaps
G = inv_lib.make_sb_matrix(ifgdates)
ixs_inc_gap = np.where(G.sum(axis=0)==0)[0]
if ixs_inc_gap.size!=0:
    with open(gapfile, 'w') as f:
        for ix in ixs_inc_gap:
            print("{}_{}".format(imdates[ix], imdates[ix+1]), file=f)


