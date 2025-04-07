import lics_unwrap as lu
import numpy as np
import os, glob

def merge_frames_pair_ext(frames, pair, outfile, ifgcoredir = 'GEOC', ext='.geo.azi.tif', fix_offset=True):
    ''' This will merge pair files of listed frames, optionally solving for their offset (not tested on phase).
    '''
    #
    '''useful hack:
    >>> tifs=glob.glob('106D*/GEOC/202503*/*azi.tif')
    >>> frames=[]
    >>> for t in tifs:
    ...     frames.append(t.split('/')[0])
    ... 
    >>> frames
    ['106D_06844_131313', '106D_07221_121306', '106D_06644_131313', '106D_07044_131313', '106D_06445_131313']
    '''
    frames = list(set(frames))
    frames.sort()
    # assuming always having overlap!
    outcore = 'tmp.mrg.'
    mediffs=[]
    tifstomerge=[]
    for i in range(len(frames)-1):
        fr1 = frames[i]
        fr2 = frames[i+1]
        fr1f = os.path.join(fr1,ifgcoredir,pair,pair+ext)
        fr2f = os.path.join(fr2, ifgcoredir, pair, pair + ext)
        b = lu.load_tif2xr(fr2f)
        b = b.where(b != 0)
        if fix_offset:
            a = lu.load_tif2xr(fr1f)
            a = a.where(a != 0)
            ba=b.interp_like(a,method='nearest')
            ba=ba.where(ba!=0)
            c=ba-a
            mediff = float(c.median())
            mediffs.append(mediff)
            b = b - np.sum(mediffs) # works ok
            # now store the corrected
            if i==0:
                outt = outcore + pair + '.' + str(i) + ext
                lu.export_xr2tif(a, outt)
                tifstomerge.append(outt)
            outt = outcore+pair+'.'+str(i+1)+ext
            lu.export_xr2tif(b,outt)
            tifstomerge.append(outt)
        else:
            if i==0:
                a = lu.load_tif2xr(fr1f)
                a = a.where(a != 0)
                outt = outcore + pair + '.' + str(i) + ext
                lu.export_xr2tif(a, outt)
                tifstomerge.append(outt)
            outt = outcore + pair + '.' + str(i + 1) + ext
            lu.export_xr2tif(b, outt)
            tifstomerge.append(outt)
    strtifs = ''
    for t in tifstomerge:
        strtifs = strtifs+t+' '
    cmd = 'gdal_merge.py -o {0} -co COMPRESS=DEFLATE -n nan -a_nodata nan {1}'.format(outfile, strtifs)
    os.system(cmd)
    print('offsets were: ')
    print(mediffs)
    print('cleaning temp files')
    os.system('rm '+strtifs)