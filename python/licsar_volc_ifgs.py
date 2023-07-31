#!/usr/bin/env python
#ML 2023
import sys, os
import volcdb

try:
    volcid=sys.argv[1]
except:
    print('Usage: licsar_volc_ifgs.py volcid [ifgs.list]')
    print('volcid ..... ID of volcano based on Powell list.')
    print('             see https://volcano.si.edu/search_volcano.cfm')
    print('             or import volcdb; volcdb.find_volcano_by_name()')
    print(' ')
    print('This will process all hires frames related to the given volcano, creating standard set of interferograms.')
    print('If you want to add your own ifgs.list with custom connections (e.g. 20210101_20210202 etc.), then ')
    print('you may check "bjobs" and when all related jobs are finished, you can rerun (in the $BATCH_CACHE_DIR/subsets/xxx/xxx) the command:')
    print('framebatch_gapfill.sh -l -P -i ifgs.list -o 5 180')
    exit()

'''
try:
    ifglist = sys.argv[2]
    print('will use file '+ifglist+' as ifg.list')
except:
    ifglist = ''
'''

volclips = volcdb.get_volclip_vids(volcid)
if not volclips:
    print('no volclips found. Check your volc id number (or contact Lin or Milan)')
    exit()
else:
    print('found {} related volclip(s).'.format(str(len(volclips))))


for volclip in volclips:
    print('processing volclip '+str(volclip))
    print('results will be stored in your $BATCH_CACHE_DIR/subsets/'+str(volclip)+' directory')
    volclipath = os.path.join(os.environ['LiCSAR_procdir'],'subsets','volc',str(volclip))
    for frameid in os.listdir(volclipath):
        print('-----')
        print('processing frame '+frameid)
        cmd = 'subset_mk_ifgs.sh '+os.path.join(volclipath, frameid) #+' '+ifglist)
        print('running command:')
        print(cmd)
        os.system(cmd)
