################################################################################
#Imports
################################################################################
import os
import sys
import numpy as np
import subprocess as subp
import scipy.spatial as spat
#import gamma functions and global configuration
from gamma_functions import *
import global_config as gc
from LiCSAR_lib.coreg_lib import get_mli_size

################################################################################
# Do unwrapping function
################################################################################
def do_unwrapping(origmasterdate,framename,ifg,ifgdir,procdir,lq,job_id):
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
    print('Filtering...')
    if not adf(difffile,filtfile,ccfile,width,'1','32',logfilename):
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
    nntree = spat.cKDTree(datapoints,leafsize=10)
    distq, gridix = nntree.query(dq)
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
    n_edges = np.histogram(abs(grid_edges),list(range(0,n_edge+1)))[0]
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
# Unwrap interferogram function
################################################################################
def unwrap_ifg(ifg, date, ifgdir, coh, width, procdir):
############################################################ Create a temp dir
    #tmpdir = os.path.join(procdir,'unwrap_tmp')
    tmpdir = os.path.join(ifgdir, date, 'unwrap_tmp')
    if os.path.exists(tmpdir) and os.path.isfile(tmpdir):
        os.remove(tmpdir)
    
    if not os.path.exists(tmpdir):
        os.mkdir(tmpdir)

############################################################ Unwrap parameters
    #parts of the interferogram with low coherence
    zeroix = coh < 0.5    
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

############################################################ Create Snaphu files
    with open(os.path.join(tmpdir,'snaphu.costinfile'),'w') as f:
        rowcost.tofile(f)
        colcost.tofile(f)

    ifg_float = np.zeros((ifg.shape[0],ifg.shape[1]*2),dtype=np.float32)
    ifg_float[:,::2] = ifg.real
    ifg_float[:,1::2] = ifg.imag
    with open(os.path.join(tmpdir,'snaphu.in'),'w') as f:
        ifg_float.tofile(f)

    with open(os.path.join(tmpdir,'snaphu.conf'),'w') as f:
        for l in get_snaphu_conf(tmpdir):
            f.write(l)

############################################################ Unwrap IFG
    logfilename = os.path.join(tmpdir,'snaphu.log')
    with open(logfilename,'w') as f:
        cmdstr = 'snaphu -v -d -f {0} {1}'.format(os.path.join(tmpdir,'snaphu.conf'),width)
        job = subp.Popen(cmdstr,shell=True,stdout=f,stderr=f,stdin=None)
    job.wait()
    try:
        uw = np.fromfile(os.path.join(tmpdir,'snaphu.out'),dtype=np.float32).reshape(ifg.shape)
    except:
        print('Something went wrong unwrapping the interferogram. Log file {0}'.format(logfilename))
        return False
    uw[coh < 0.5] = np.nan
    try:
        uw.byteswap().tofile(os.path.join(ifgdir,date,date+'.unw'))
    except:
        print('Something went wrong writing the unwrapped interferogram to file {0}.'.format(os.path.join(ifgdir,date+'.unw')))
        return False


    return True
    #cleaning the tmpdir
    if os.path.exists(tmpdir) and os.path.isfile(tmpdir):
        os.remove(tmpdir)

                             

################################################################################
# Get snaphu conf
################################################################################
def get_snaphu_conf(tmpdir):
    return ('INFILE  {0}/snaphu.in\n'.format(tmpdir),
            'OUTFILE {0}/snaphu.out\n'.format(tmpdir),
            'COSTINFILE {0}/snaphu.costinfile\n'.format(tmpdir),
            'STATCOSTMODE  DEFO\n',
            'INFILEFORMAT  COMPLEX_DATA\n',
            'OUTFILEFORMAT FLOAT_DATA\n',
            'NSHORTCYCLE 100\n')
