#!/usr/bin/env python3
# two outputs are: png plot and text file for gaps (both are mandatory)
# e.g. $LiCSAR_public/$track/$frame output.png gaps.txt
#%% Import
import os
import sys
import numpy as np
import LiCSBAS_io_lib as io_lib
import LiCSBAS_plot_lib as plot_lib
import LiCSBAS_tools_lib as tools_lib
import LiCSBAS_inv_lib as inv_lib

#%% File setting
framedir = sys.argv[1]
pngfile = sys.argv[2] 
gapfile = sys.argv[3]

ifgdir = os.path.join(framedir, 'products')
bperp_file = os.path.join(framedir, 'metadata', 'baselines')

#%%
ifgdates = tools_lib.get_ifgdates(ifgdir)
imdates = tools_lib.ifgdates2imdates(ifgdates)

if not os.path.exists(bperp_file):
    print('Make dummy bperp')
    bperp_file = 'baselines_tmp.txt'
    io_lib.make_dummy_bperp(bperp_file, imdates)
elif not io_lib.read_bperp_file(bperp_file, imdates):
    print('Make dummy bperp')
    bperp_file = 'baselines_tmp.txt'
    io_lib.make_dummy_bperp(bperp_file, imdates)

bperp = io_lib.read_bperp_file(bperp_file, imdates)

plot_lib.plot_network(ifgdates, bperp, [], pngfile)

## Identify gaps
G = inv_lib.make_sb_matrix(ifgdates)
ixs_inc_gap = np.where(G.sum(axis=0)==0)[0]
if ixs_inc_gap.size!=0:
    with open(gapfile, 'w') as f:
        for ix in ixs_inc_gap:
            print("{}_{}".format(imdates[ix], imdates[ix+1]), file=f)

