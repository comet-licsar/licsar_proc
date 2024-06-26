################################################################################
#Imports
################################################################################
import os
import shutil
import sys
import re
import numpy as np
import subprocess as subp
import scipy.spatial as spat
#import gamma functions and global configuration
from gamma_functions import *
import global_config as gc
from LiCSAR_lib.coreg_lib import get_mli_size

################################################################################
# Do unwrapping function (lq is needed only if job_id is set...)
################################################################################
def do_unwrapping(origmasterdate,ifg,ifgdir,procdir,lq,job_id):
############################################################ Get master mli, length and width
    mastermli = os.path.join(procdir,'SLC',
                             origmasterdate,
                             origmasterdate+
                             '.slc.mli')
                             
    #read width and length
    [width, length] = get_mli_size(mastermli+'.par')

    print('\nUnwrapping combination {0}'.format(ifg))

############################################################ Filter using an adaptive filter algorithm 
    ifgdirthis = os.path.join(ifgdir,ifg)
    difffile = os.path.join(ifgdirthis,ifg+'.diff')
    filtfile = os.path.join(ifgdirthis,ifg+'.filt.diff')
    ccfile = os.path.join(ifgdirthis,ifg+'.filt.cc')
    logfilename = os.path.join(procdir,'log','adf_{0}.log'.format(ifg))
    if os.path.exists(filtfile):
        print('filtered ifg exists - will use it further')
    else:
        print('Filtering...')
        if not adf(difffile,filtfile,ccfile,width,str(gc.adf_alpha),str(gc.adf_window),logfilename):
            print('\nERROR:', file=sys.stderr)
            print('\nSomething went wrong filtering interferogram{0}.'.format(ifg), file=sys.stderr)
            return 1
        #Create a preview ras file
        logfilename = os.path.join(procdir,'log','rasmph_pwr_filt_{0}.log'.format(ifg))
        if not rasmph_pwr(filtfile,mastermli,width,logfilename):
            print('\nERROR:', file=sys.stderr)
            print('\nSomething went wrong creating the preview raster file of the filtered interferogram {0}.'.format(ifg), file=sys.stderr)
            return 1

############################################################ Unwrap the interferogram
    # ccfile = os.path.join(ifgdirthis,filtfile[:-5]+'.cc')
    coh = np.fromfile(ccfile,dtype=np.float32).byteswap().reshape((int(length),int(width)))
    ccthres = gc.coh_unwrap_threshold
    numoverthres = len(coh[coh>ccthres])
    numall = len(coh[coh>0])
    #now let's choose proper ratio to select unwrapper
    if numall == 0 or numoverthres == 0:
        print('\n no points above coh threshold to unwrap! cancelling as error', file=sys.stderr)
        return 1
    elif numall/numoverthres > 2.2:
        print('the number of selected points is very low. skipping internal routines and using gamma mcf')
        if not unwrap_ifg_gamma(ifg, ifgdir, width, length, coh):
            print('\nERROR:', file=sys.stderr)
            print('\nSomething went wrong unwrapping ifg {0} using gamma mcf'.format(ifg), file=sys.stderr)
            return 2
    else:
        ifg_w = np.fromfile(filtfile,dtype=np.complex64).byteswap().reshape((int(length),int(width)))
        print('Unwrapping...')
        if not unwrap_ifg(ifg_w, ifg, ifgdir, coh,width,procdir):
            print('\nERROR:', file=sys.stderr)
            print('\nSomething went wrong unwrapping interferogram {0}.'.format(ifg), file=sys.stderr)
            return 2
    
    #Create a ras file for the unwrapped interferogram
    logfilename = os.path.join(procdir,'log','rasrmg_{0}.log'.format(ifg))
    if not rasrmg(os.path.join(ifgdir,ifg,ifg+'.unw'),mastermli,width,1,logfilename):
        print('\nERROR:', file=sys.stderr)
        print('\nSomething went wrong creating the preview raster file of the unwrapped interferogram {0}.'.format(ifg), file=sys.stderr)
        return 1

    # Log the update of the IFG to indicate it is now unwrapped in the DB
    if job_id != -1:
        lq.update_ifg_product_unwrapped(job_id, ifg+'.diff')

    return 0
 
################################################################################
# Get edges function
################################################################################
def get_edges(ph,zeroix):
    length, width = ph.shape
    i,j = np.where(~zeroix)
    datapoints = np.array((i,j)).T
    ix = np.ones(ph.shape,dtype=bool)
    iq,jq = np.where(ix)
    dq = np.array((iq,jq)).T
    nntree = spat.cKDTree(datapoints,leafsize=10,compact_nodes=False,balanced_tree=False)
    distq, gridix = nntree.query(dq,n_jobs=-1)
    rowedges = np.array((gridix[:-width],
                         gridix[width:])).T
    gridixT = np.reshape(gridix,(length,width)).T.flatten()
    coledges = np.array((gridixT[:-length],
                         gridixT[length:])).T
    grid_edges = np.vstack((rowedges,coledges))
    sortix = np.argsort(grid_edges, axis=1)
    sort_edges = np.sort(grid_edges,axis=1)
    edge_sign = sortix[:,1]-sortix[:,0]
    sort_edges_flat = sort_edges[:,0]+sort_edges[:,1]*sort_edges.shape[0]
    dummy, alledge_ix, inverse_ix = np.unique(sort_edges_flat,
                                              return_index=True, 
                                              return_inverse=True)
    alledges = sort_edges[alledge_ix,:]
    alledges[alledges[:,0] == alledges[:,1]] = -1
    alledges_flat = alledges[:,0]+alledges[:,1]*alledges.shape[0]
    dummy, edge_ix, inverse_ix2 = np.unique(alledges_flat,
                                            return_index=True, 
                                            return_inverse=True)
    edges = alledges[edge_ix,:]
    edges = np.delete(edges,0,axis=0)
    n_edge = edges.shape[0]
    edges = np.hstack((np.arange(n_edge)[:,None]+1,edges))
    gridedgeix = (inverse_ix2[inverse_ix]-1)
    sameixval = gridedgeix.max()+1
    gridedgeix[gridedgeix == -1] = sameixval
    gridedgeix *= edge_sign
    rowix = np.ma.array(np.reshape(gridedgeix[:width*
                                              (length-1)],
                       (length-1,width)))
    rowix = np.ma.masked_where(abs(rowix) == sameixval,rowix)
    colix = np.ma.array(np.reshape(gridedgeix[width*
                                              (length-1):],
                       (width-1,length)).T)
    colix = np.ma.masked_where(abs(colix) == sameixval,colix)
    return edges, n_edge, rowix, colix, gridix

################################################################################
# Get costs
################################################################################
def get_costs(edges, n_edge, rowix, colix, zeroix):
    length, width = zeroix.shape
    maxshort = 32000
    costscale = 40
    nshortcycle = 100
    i,j = np.where(~zeroix)
    grid_edges = np.vstack((rowix.compressed()[:,None],colix.compressed()[:,None]))
    n_edges =  np.histogram(abs(grid_edges),n_edge,(0,n_edge))[0]
    edge_length = np.sqrt(np.diff(i[edges[:,1:]],axis=1)**2+
                          np.diff(j[edges[:,1:]],axis=1)**2)
    sigsq_noise = np.zeros(edge_length.shape,dtype=np.float32)
    aps_range = 20000
    sigsq_aps = (2*np.pi)**2
    sigsq_noise += sigsq_aps*(1-np.exp(-edge_length
                                       * 80
                                       * 3 / aps_range))
    sigsq_noise /= 10
    sigsq = np.int16((sigsq_noise*nshortcycle**2)/costscale*n_edges[:,None])
    sigsq[sigsq<1] = 1
    rowcost = np.zeros((length-1,width*4),dtype=np.int16)
    colcost = np.zeros((length,(width-1)*4),dtype=np.int16)
    rowstdgrid = np.ones(rowix.shape,dtype=np.int16)
    colstdgrid = np.ones(colix.shape,dtype=np.int16)
    rowcost[:,2::4] = maxshort
    colcost[:,2::4] = maxshort
    rowcost[:,3::4] = (-1-maxshort)+1
    colcost[:,3::4] = (-1-maxshort)+1
    if hasattr(rowix.mask,"__len__"):
        rowstdgrid[~rowix.mask] = sigsq[abs(rowix.compressed())].squeeze()
    else:
        mask = np.ones(rowstdgrid.shape,dtype=np.bool_)
        rowstdgrid[mask] = sigsq[abs(rowix.compressed())].squeeze()
    rowcost[:,1::4] = rowstdgrid
    if hasattr(colix.mask,"__len__"):
        colstdgrid[~colix.mask] = sigsq[abs(colix.compressed())].squeeze()
    else:
        mask = np.ones(colstdgrid.shape,dtype=np.bool_)
        colstdgrid[mask] = sigsq[abs(colix.compressed())].squeeze()
    colcost[:,1::4] = colstdgrid
    return rowcost, colcost

################################################################################
# Unwrap interferogram function(s)
################################################################################
def unwrap_ifg_gamma(date, ifgdir, width, length, coh):
    #use mcf with coh and mask above threshold
    ifgfile = os.path.join(ifgdir,date,date+'.filt.diff')
    ccfile = os.path.join(ifgdir,date,date+'.filt.cc')
    rc = os.system('rascc_mask {0} - {1} 1 1 0 1 1 {2} 0.1 0.9 1 0.35 1 *.bmp'.format(ccfile, str(width), str(gc.coh_unwrap_threshold)))
    maskfile = ccfile+'_mask.bmp'
    if not os.path.exists(maskfile):
        print('ERROR during generating mask file - maybe not saved as BMP?')
        return False
    unwfile = os.path.join(ifgdir,date,date+'.unw')
    if os.path.exists(unwfile):
        print('ERROR: actually the unwrapped file seems to exist! not overwriting, cancelling')
        return False
    rc = os.system('mcf {0} {1} {2} {3} {4}'.format(ifgfile, ccfile, maskfile, unwfile, str(width)))
    if not os.path.exists(unwfile):
        print('ERROR in unwrapping')
        return False
    #now just to (really) mask the unw file..
    unw2 = np.fromfile(unwfile,dtype=np.float32).byteswap().reshape((int(length),int(width)))
    #maybe not needed, but keeping nans, as in snaphu approach
    unw[coh <= gc.coh_unwrap_threshold] = np.nan
    #unw = unw - np.nanmedian(unw)  #actually should be already
    try:
        os.remove(unwfile)
        unw.byteswap().tofile(unwfile)
    except:
        print('some error saving masked unw file')
        return False
    return True


def unwrap_ifg(ifg, date, ifgdir, coh, width, procdir):
    #original function, to use snaphu..
############################################################ Create a temp dir
    #tmpdir = os.path.join(procdir,'unwrap_tmp')
    tmpdir = os.path.join(ifgdir, date, 'unwrap_tmp')
    if os.path.exists(tmpdir) and os.path.isfile(tmpdir):
        os.remove(tmpdir)
    if not os.path.exists(tmpdir):
        os.mkdir(tmpdir)
    maskfile = os.path.join(tmpdir,'snaphu.mask')
############################################################ Unwrap parameters
    #parts of the interferogram with low coherence
    zeroix = coh < gc.coh_unwrap_threshold  #orig value was 0.5
    #zero low coherence parts of the interferogram
    ifg[zeroix] = np.complex64(0+0j)
    #get edges
    edges, n_edge, rowix, colix, gridix = get_edges(ifg, zeroix)
    #create costs
    rowcost, colcost = get_costs(edges, n_edge, rowix, colix, zeroix)
    #find spots which are not below coherence threshold
    ithis, jthis = np.where(~zeroix)
    #patch zeroed parts of the IFG?
    ifg = ifg[ithis[gridix],jthis[gridix]].reshape(ifg.shape)
    #prepare mask
    mask = (coh > gc.coh_unwrap_threshold).astype(np.byte)
    mask.tofile(maskfile)
############################################################ Create Snaphu files
    with open(os.path.join(tmpdir,'snaphu.costinfile'),'w') as f:
        rowcost.tofile(f)
        colcost.tofile(f)
    ifg_float = np.zeros((ifg.shape[0],ifg.shape[1]*2),dtype=np.float32)
    ifg_float[:,::2] = ifg.real
    ifg_float[:,1::2] = ifg.imag
    with open(os.path.join(tmpdir,'snaphu.in'),'w') as f:
        ifg_float.tofile(f)
    #copy customized snaphu.conf if exists for the given frame
    if os.path.exists(os.path.join(procdir,'log','snaphu.conf')):
        shutil.copyfile(os.path.join(procdir,'log','snaphu.conf'), os.path.join(tmpdir,'snaphu.conf'))
        subp.call(["sed", "-i", 's/TMPDIR/{}/'.format(re.sub('/','\/',tmpdir)), os.path.join(tmpdir,'snaphu.conf')])
    else:
        with open(os.path.join(tmpdir,'snaphu.conf'),'w') as f:
            for l in get_snaphu_conf(tmpdir):
                f.write(l)
############################################################ Unwrap IFG
    logfilename = os.path.join(tmpdir,'snaphu.log')
    with open(logfilename,'w') as f:
        #new snaphu: parameter -S for optimization within tiling
        #cmdstr = 'snaphu -v -S -d -f {0} {1}'.format(os.path.join(tmpdir,'snaphu.conf'),width)
        #old snaphu (YuM identified significant differences, the older version is preferred)
        #cmdstr = 'snaphu -v -d -f {0} {1}'.format(os.path.join(tmpdir,'snaphu.conf'),width)
        #including coh-based mask
        cmdstr = 'snaphu -v -d -M {0} -f {1} {2}'.format(maskfile, os.path.join(tmpdir,'snaphu.conf'),width)
        job = subp.Popen(cmdstr,shell=True,stdout=f,stderr=f,stdin=None)
    job.wait()
    try:
        uw = np.fromfile(os.path.join(tmpdir,'snaphu.out'),dtype=np.float32).reshape(ifg.shape)
    except:
        print('Something went wrong unwrapping the interferogram. Log file {0}'.format(logfilename))
        return False
    # fully removing (filtered) coh points of < coh threshold
    # hmm... keeping it customisable..
    uw[coh < gc.coh_unwrap_threshold] = np.nan
    try:
        uw.byteswap().tofile(os.path.join(ifgdir,date,date+'.unw'))
    except:
        print('Something went wrong writing the unwrapped interferogram to file {0}.'.format(os.path.join(ifgdir,date+'.unw')))
        return False
    #cleaning the tmpdir
    #if os.path.exists(tmpdir) and os.path.isfile(tmpdir):
    #    os.remove(tmpdir)
    if os.path.exists(tmpdir) and os.path.isdir(tmpdir):
        shutil.rmtree(tmpdir)
    return True

def unwrap_geo(procdir, frame, pair):
    rc = os.system('cd {0}; unwrap_geo.sh {1} {2}'.format(procdir, frame, pair))
    if os.path.exists(os.path.join(procdir,'GEOC',pair,pair+'.geo.unw.tif')):
        rc = 0
    else:
        rc = 1
    return rc


################################################################################
# Get snaphu conf
################################################################################
# constant problems with tiling - returning back to no-tiling
def get_snaphu_conf(tmpdir):
    return ('INFILE  {0}/snaphu.in\n'.format(tmpdir),
            'OUTFILE {0}/snaphu.out\n'.format(tmpdir),
            'COSTINFILE {0}/snaphu.costinfile\n'.format(tmpdir),
            'STATCOSTMODE  DEFO\n',
            'INFILEFORMAT  COMPLEX_DATA\n',
            'OUTFILEFORMAT FLOAT_DATA\n',
            'NSHORTCYCLE 100\n',
            'NTILEROW 1\n',
            'NTILECOL 1\n',
            'NPROC 1\n',
            'RMTMPTILE TRUE\n')

def get_snaphu_conf_tiled(tmpdir):
    return ('INFILE  {0}/snaphu.in\n'.format(tmpdir),
            'OUTFILE {0}/snaphu.out\n'.format(tmpdir),
            'COSTINFILE {0}/snaphu.costinfile\n'.format(tmpdir),
            'STATCOSTMODE  DEFO\n',
            'INFILEFORMAT  COMPLEX_DATA\n',
            'OUTFILEFORMAT FLOAT_DATA\n',
            'NSHORTCYCLE 100\n',
            'NTILEROW 4\n',
            'NTILECOL 4\n',
            'ROWOVRLP 200\n',
            'COLOVRLP 200\n',
            'NPROC 16\n',
            'RMTMPTILE TRUE\n')


def demedian_unw(filename):
    #this is to reduce the unw by a median
    a = np.fromfile(filename,dtype=np.float32)
    print('reducing by median')
    a = a - np.nanmedian(a)
    os.remove(filename)
    a.tofile(filename)



def reunwrap(frame, pair, procdir = os.getcwd()):
    '''Simple function to call the updated unwrapper and regenerated the unw, by rewriting it..

    e.g. frame='106D_05447_131313', pair='20210414_20210520'
    '''
    track = str(int(frame[:3]))
    ml = 1
    from lics_unwrap import process_ifg, export_xr2tif, wrap2phase
    ifg = process_ifg(frame, pair, procdir, ml=ml, fillby='nearest', lowpass=False, specmag=True, )
    outdir=os.path.join(os.environ['LiCSAR_public'], track, frame,'interferograms', pair)
    outunw = os.path.join(outdir, pair+'.geo.unw.tif')
    outunwpng = os.path.join(outdir, pair+'.geo.unw.png')
    outpha = os.path.join(outdir, pair+'.geo.diff_pha.tif')
    outphapng = os.path.join(outdir, pair+'.geo.diff.png')
    if os.path.exists(outphapng):
        os.remove(outphapng)
        os.remove(outunwpng)
    #
    os.remove(outpha)
    os.remove(outunw)
    ifg['origpha'].values = wrap2phase(ifg.toremove+ifg.filtpha)
    export_xr2tif(ifg['origpha'], outpha)
    export_xr2tif(ifg['unw'], outunw)
    #
    # creating previews
    cmd = "x="+pair+"; source $LiCSARpath/lib/*; create_preview_wrapped "+outdir+"/$x/$x.geo.diff_pha.tif"
    rc = os.system(cmd)
    cmd = "x="+pair+"; source $LiCSARpath/lib/*; create_preview_unwrapped "+outdir+"/$x/$x.geo.unw.tif"
    rc = os.system(cmd)
