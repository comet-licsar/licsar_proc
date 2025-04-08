import lics_unwrap as lu
import numpy as np
import os, glob
import framecare as fc

def merge_frames_pair_ext(frames, pair = None, outfile = 'merged.tif', ifgcoredir = 'GEOC', ext='.geo.azi.tif', fix_offset=True,
                          mask_by_ref = False, usepublic = False):
    ''' This will merge files of listed frames, optionally solving for their offset (not tested on phase).
    If pair is given, will find data under pair subfolder, otherwise it tries find directly in ifgcoredir.
    in case of 'usepublic', the ifgcoredir can be None, or one of pubdirs such as 'interferograms', 'metadata'..
    '''
    #
    '''useful hack:
    >>> tifs=glob.glob('106D*/GEOC/202503*/*azi.tif')
frames=[]
for t in tifs:
    frames.append(t.split('/')[0])
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
    if not pair:
        ispair = False
        pair = ''
    else:
        ispair = True
    if not ifgcoredir:
        usepublic = True
        if ispair:
            ifgcoredir = 'interferograms'
        else:
            ifgcoredir = 'metadata'
    tr=str(int(frames[0][:3]))
    for i in range(len(frames)-1):
        fr1 = frames[i]
        fr2 = frames[i+1]
        if ispair:
            if usepublic:
                fr1dir = os.path.join(os.environ['LiCSAR_public'], tr, fr1)
                fr2dir = os.path.join(os.environ['LiCSAR_public'], tr, fr2)
                fr1f = os.path.join(fr1dir, ifgcoredir, pair, pair + ext)
                fr2f = os.path.join(fr2dir, ifgcoredir, pair, pair + ext)
            else:
                fr1f = os.path.join(fr1, ifgcoredir, pair, pair + ext)
                fr2f = os.path.join(fr2, ifgcoredir, pair, pair + ext)
        else:
            if usepublic:
                fr1dir = os.path.join(os.environ['LiCSAR_public'], tr, fr1)
                fr2dir = os.path.join(os.environ['LiCSAR_public'], tr, fr2)
                fr1f = os.path.join(fr1dir, ifgcoredir, fr1 + ext)
                fr2f = os.path.join(fr2dir, ifgcoredir, fr2 + ext)
            else:
                fr1f = os.path.join(fr1,ifgcoredir,fr1+ext)
                fr2f = os.path.join(fr2, ifgcoredir, fr2+ ext)
        if mask_by_ref:
            m1 = fc.get_master(fr1)
            reff1 = os.path.join(os.environ['LiCSAR_public'], tr, fr1, 'epochs', m1, m1+'.geo.mli.tif')
            ref1 = lu.load_tif2xr(reff1)
            ref1 = ref1.where(ref1 != 0)
            m2 = fc.get_master(fr2)
            reff2 = os.path.join(os.environ['LiCSAR_public'], tr, fr2, 'epochs', m2, m2+'.geo.mli.tif')
            ref2 = lu.load_tif2xr(reff2)
            ref2 = ref2.where(ref2 != 0)
        b = lu.load_tif2xr(fr2f)
        if mask_by_ref:
            b = b.where(~np.isnan(ref2))
        b = b.where(b != 0)
        if fix_offset:
            a = lu.load_tif2xr(fr1f)
            if mask_by_ref:
                a = a.where(~np.isnan(ref1))
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
                if mask_by_ref:
                    a = a.where(~np.isnan(ref1))
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


def merge_enus_frames(frames, outcore = None):
    ''' This would merge enus for frames in list, saving as [outcore].geo.E.tif etc.
    -- if outcore == None, it will use the frame orb name '''
    if type(outcore) == type(None):
        outcore = frames[0].split('_')[0]
    for enu in ['E', 'N', 'U', 'hgt']:
        ext='.geo.{0}.tif'.format(enu)
        outfile = outcore+ext
        merge_frames_pair_ext(frames, pair=None, outfile=outfile, ifgcoredir='metadata', ext=ext,
                              fix_offset=False, mask_by_ref=True, usepublic=True)


