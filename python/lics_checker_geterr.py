#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Nov  5 09:34:59 2019

@author: earymo
"""

#%%
import pickle
import sys
import os
import shutil
import numpy as np

#%%
licscheckdir = '/nfs/a1/insar/lics_check/'
if not os.path.exists(licscheckdir):
    licscheckdir = '.'

if len(sys.argv) <= 1:
    print('\nUsage: {} FrameID'.format(os.path.basename(sys.argv[0])))
    print('Example: {} 138D_05321_080913\n'.format(os.path.basename(sys.argv[0])))
    print('{}[FrameID]errtmp which contains png with errors will be created.\n'.format(licscheckdir))
    sys.exit(1)


#%%
frame = sys.argv[1]

results = os.path.join(licscheckdir, '{}.savedResults'.format(frame))
if not os.path.exists(results):
    print('No {} exists!'.format(results))
    sys.exit(1)

errdir = os.path.join(licscheckdir, '{}errtmp'.format(frame))

if os.path.exists(errdir):
    shutil.rmtree(errdir)
os.mkdir(errdir)


#%%
data = pickle.load(open(results, 'rb'))
dates = np.array(data[0])
errors = np.array(data[1])

ixs_err = np.where(errors!=0)[0]

print('Error codes:\n1 = swath discontinuity\n2 = visible tiles (unwrapping error)\n3= coregistration error (ghosts, shifted areas)\n4= missing bursts (grey/blue areas)\n5= missing unwrapped data\n6= other error\n7= ESD error\n')
for ix_err in ixs_err:
    print('{}  {}'.format(dates[ix_err], errors[ix_err]))
    unwpng = dates[ix_err]+'.geo.unw.png'
    diffpng = dates[ix_err]+'.geo.diff.png'
    unwpng_src = os.path.join(licscheckdir, frame, unwpng)
    diffpng_src = os.path.join(licscheckdir, frame, diffpng)
    unwpng_dst = os.path.join(errdir, unwpng)
    diffpng_dst = os.path.join(errdir, diffpng)
    os.symlink(unwpng_src, unwpng_dst)
    os.symlink(diffpng_src, diffpng_dst)

print('{}/{} ifgs with errors'.format(len(ixs_err), len(errors)))
print('\nTo view png with errors:')
print('eog {} &\n'.format(errdir))
