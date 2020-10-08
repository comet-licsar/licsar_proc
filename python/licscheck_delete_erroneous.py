#!/usr/bin/env python3

import pickle
import sys
import os
import shutil
import numpy as np

if len(sys.argv) <= 1:
    print('\nUsage: {} savedResults_file'.format(os.path.basename(sys.argv[0])))
    print('Example: {} 138D_05321_080913.savedResults\n'.format(os.path.basename(sys.argv[0])))
    print('errtmp which contains png with errors will be created.\n')
    sys.exit(1)


#%%
results = sys.argv[1]
frame = results.split('.')[0]

if not os.path.exists(results):
    print('No {} exists!'.format(results))
    sys.exit(1)

data = pickle.load(open(results, 'rb'))
dates = np.array(data[0])
errors = np.array(data[1])
ixs_err = np.where(errors!=0)[0]

print('Error codes:\n1 = swath discontinuity\n2 = visible tiles (unwrapping error)\n3= coregistration error (ghosts, shifted areas)\n '\
      '4= missing bursts (grey/blue areas)\n5= missing unwrapped data\n6= other error\n7= ESD error\n')
for ix_err in ixs_err:
    print('{}  {}'.format(dates[ix_err], errors[ix_err]))
    #unwpng = dates[ix_err]+'.geo.unw.png'
    #diffpng = dates[ix_err]+'.geo.diff.png'
    #unwpng_src = os.path.join(frame, unwpng)
    #diffpng_src = os.path.join(frame, diffpng)
    #unwpng_dst = os.path.join(errdir, unwpng)
    #diffpng_dst = os.path.join(errdir, diffpng)
    #os.symlink(unwpng_src, unwpng_dst)
    #os.symlink(diffpng_src, diffpng_dst)

print('{}/{} ifgs with errors'.format(len(ixs_err), len(errors)))

errors_to_delete = [1, 2, 3, 7]
errors_to_delete = [7]
epochs = set()
print('removing those of errors {}'.format(errors_to_delete))
for error in errors_to_delete:
    ixs_err = np.where(errors==error)[0]
    for ix_err in ixs_err:
        print(dates[ix_err])
        for epoch in dates[ix_err].split('_'):
            epochs.add(epoch)
        os.system('remove_from_lics.sh {0} {1}'.format(frame,dates[ix_err]))

procdir = '/gws/nopw/j04/nceo_geohazards_vol1/projects/LiCS/proc/current'
tr=str(int(frame[:3]))
for epoch in epochs:
    lut = os.path.join(procdir,tr,frame,'LUT',epoch+'.7z')
    rslc = os.path.join(procdir,tr,frame,'RSLC',epoch+'.7z')
    for z in [lut,rslc]:
        if os.path.exists(z):
            os.remove(z)
