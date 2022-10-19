#!/usr/bin/env python3
# this is to decompose two unw ifgs - very basic way
import subprocess as subp
import numpy as np
from scipy import interpolate
import interseis_lib as lib
import matplotlib.pyplot as plt
import importlib

# these packages are only needed for the final multivariate plot
import seaborn as sns
import pandas as pd

importlib.reload(lib)



# 2022-10-18 starts here

frame_desc = '082D_05128_030500'
frame_asc = '002A_05136_020502'
asctif='saojorge_offset_A.tif'
desctif='saojorge_offset_D.tif'
beta=25.83

inc_asc, heading_asc = get_frame_inc_heading(frame_asc)
inc_desc, heading_desc = get_frame_inc_heading(frame_desc)
asc = load_tif2xr(asctif)
desc = load_tif2xr(desctif)

# now transform it towards asc:
cube=xr.Dataset()
cube['asc'] = asc
cube['desc'] = desc.interp_like(asc, method='linear'); desc=None
cube['asc_inc'] = inc_asc.interp_like(asc, method='linear'); inc_asc=None
cube['desc_inc'] = inc_desc.interp_like(asc, method='linear'); inc_desc=None
cube['asc_heading'] = heading_asc.interp_like(asc, method='linear'); heading_asc=None
cube['desc_heading'] = heading_desc.interp_like(asc, method='linear'); heading_desc=None

# and decompose:
cube['U']=cube.asc.copy()
cube['E']=cube.asc.copy()
cube['U'].values, cube['E'].values = decompose_np(cube.asc, cube.desc, cube.asc_heading, cube.desc_heading, cube.asc_inc, cube.desc_inc)

from lics_unwrap import *
def get_frame_inc_heading(frame):
    geoframedir = os.path.join(os.environ['LiCSAR_public'], str(int(frame[:3])), frame)
    # look angle (inc) / heading - probably ok, but needs check:
    e=os.path.join(geoframedir,'metadata',frame+'.geo.E.tif')
    #n=os.path.join(geoframedir,'metadata',frame+'.geo.N.tif') #no need for N
    u=os.path.join(geoframedir,'metadata',frame+'.geo.U.tif')
    e = load_tif2xr(e)
    #n = load_tif2xr(n, cliparea_geo=cliparea)
    u = load_tif2xr(u)
    #
    theta=np.arcsin(u)
    phi=np.arccos(e/np.cos(theta))
    heading = np.rad2deg(phi)-180
    inc = 90-np.rad2deg(theta)   #correct
    #inc.values.tofile(outinc)
    return inc, heading


import dask.array as da
def decompose_dask(cube, blocklen=5, num_workers=5):
    winsize = (blocklen, blocklen)
    asc = da.from_array(cube['asc'].astype(np.float32), chunks=winsize)
    desc = da.from_array(cube['desc'].astype(np.float32), chunks=winsize)
    ascinc = da.from_array(cube['asc_inc'].astype(np.float32), chunks=winsize)
    descinc = da.from_array(cube['desc_inc'].astype(np.float32), chunks=winsize)
    aschead = da.from_array(cube['asc_heading'].astype(np.float32), chunks=winsize)
    deschead = da.from_array(cube['desc_heading'].astype(np.float32), chunks=winsize)
    #f = da.map_blocks(decompose_np, asc, desc, aschead, deschead, ascinc, descinc, beta=0, meta=np.array((),())) #, chunks = (1,1))
    f = da.map_blocks(decompose_np, asc, desc, aschead, deschead, ascinc, descinc, beta=0, meta=(np.array((), dtype=np.float32), np.array((), dtype=np.float32)))
    return f.compute(num_workers=num_workers)

# 2022-10-18 - this should be pretty good one (next only use weights or something)
def decompose_np(vel_asc, vel_desc, aschead, deschead, ascinc, descinc, beta=0):
    '''Decomposes values from ascending and descending geotiffs, using heading and inc. angle
    (these might be arrays of same size of just float values)
    
    Args:
        beta (float): angle of expected horizontal motion direction, clockwise from the E, in degrees
    '''
    vel_E = np.zeros(vel_desc.shape)
    vel_U = np.zeros(vel_desc.shape)
    #
    U_asc = np.cos(np.radians(ascinc))
    E_asc = -np.sin(np.radians(ascinc))*np.cos(np.radians(aschead+beta))
    U_desc = np.cos(np.radians(descinc))
    E_desc = -np.sin(np.radians(descinc))*np.cos(np.radians(deschead+beta))
    #
    for ii in np.arange(0,vel_E.shape[0]):
        for jj in np.arange(0,vel_E.shape[1]):
            # velocities
            d = np.array([[vel_asc[ii,jj], vel_desc[ii,jj]]]).T
            # if the velocities contain nan, will return nan:
            if np.isnan(np.max(d)):
                vel_U[ii,jj] = np.nan
                vel_E[ii,jj] = np.nan
            else:
                # create the design matrix
                if np.isscalar(U_asc):  # in case of only values (i.e. one inc and heading per each frame)
                    G = np.array([[U_asc, E_asc], [U_desc, E_desc]])
                else:  # in case this is array
                    G = np.array([[U_asc[ii,jj], E_asc[ii,jj]], [U_desc[ii,jj], E_desc[ii,jj]]])
                # solve the linear system for the Up and East velocities
                m = np.linalg.solve(G, d)
                # save to arrays
                vel_U[ii,jj] = m[0]
                vel_E[ii,jj] = m[1]
    return vel_U, vel_E


'''

aschead=-9.918319
deschead=-169.61931
ascinc=39.4824
descinc=33.7491
desc=xr.open_dataset('082D_05128_030500ok.nc')
asc=xr.open_dataset('002A_05136_020502.nc')



D=desc.cum[-2]-desc.cum[-3]
A=asc.cum[-1]-asc.cum[-4]

from unwrp_multiscale import *
export_xr2tif(A,'A.tif', dogdal=False)
export_xr2tif(D,'D.tif', dogdal=False)
os.system('gdalwarp2match.py A.tif D.tif Aok.tif')
os.system('gdalwarp2match.py D.tif Aok.tif Dok.tif')
'''
#def decompose_xr(asc, desc, aschead, deschead, ascinc, descinc):
#    
U,E = decompose('Aok.tif', 'Dok.tif', aschead, deschead, ascinc, descinc)
aa = rioxarray.open_rasterio('Aok.tif')
aa.values[0]=U
export_xr2tif(aa,'U.tif', lonlat=False,dogdal=False)

import LiCSBAS_io_lib as io
# to load to np
asctif=...
desctif=...
vel_asc = io.read_geotiff(asctif)
vel_desc = io.read_geotiff(desctif)



aschead=349.613
deschead=190.3898
ascinc=34.403
descinc=34.240








# setup file names
vel_file_asc = 'data/087A_04904_121313_vel'
par_file_asc = 'data/087A_04904_121313.par'
E_file_asc = 'data/087A_04904_121313_E.geo'
N_file_asc = 'data/087A_04904_121313_N.geo'
U_file_asc = 'data/087A_04904_121313_U.geo'

vel_file_desc = 'data/167D_04884_131212_vel'
par_file_desc = 'data/167D_04884_131212.par'
E_file_desc = 'data/167D_04884_131212_E.geo'
N_file_desc = 'data/167D_04884_131212_N.geo'
U_file_desc = 'data/167D_04884_131212_U.geo'

# read array dimensions from par file
width_asc = int(lib.get_par(par_file_asc,'width'))
length_asc = int(lib.get_par(par_file_asc,'nlines'))

width_desc = int(lib.get_par(par_file_desc,'width'))
length_desc = int(lib.get_par(par_file_desc,'nlines'))

# get corner positions
corner_lat_asc = float(lib.get_par(par_file_asc,'corner_lat'))
corner_lon_asc = float(lib.get_par(par_file_asc,'corner_lon'))

corner_lat_desc = float(lib.get_par(par_file_desc,'corner_lat'))
corner_lon_desc = float(lib.get_par(par_file_desc,'corner_lon'))

# get post spacing (distance between velocity measurements)
post_lat_asc = float(lib.get_par(par_file_asc,'post_lat'))
post_lon_asc = float(lib.get_par(par_file_asc,'post_lon'))

post_lat_desc = float(lib.get_par(par_file_desc,'post_lat'))
post_lon_desc = float(lib.get_par(par_file_desc,'post_lon'))

# calculate grid spacings
lat_asc = corner_lat_asc + post_lat_asc*np.arange(1,length_asc+1) - post_lat_asc/2
lon_asc = corner_lon_asc + post_lon_asc*np.arange(1,width_asc+1) - post_lon_asc/2

lat_desc = corner_lat_desc + post_lat_desc*np.arange(1,length_desc+1) - post_lat_desc/2
lon_desc = corner_lon_desc + post_lon_desc*np.arange(1,width_desc+1) - post_lon_desc/2

# load in velocities
vel_asc = np.fromfile(vel_file_asc, dtype='float32').reshape((length_asc, width_asc))
vel_desc = np.fromfile(vel_file_desc, dtype='float32').reshape((length_desc, width_desc))

# load in unit vectors
E_asc = np.fromfile(E_file_asc, dtype='float32').reshape((length_asc, width_asc))
N_asc = np.fromfile(N_file_asc, dtype='float32').reshape((length_asc, width_asc))
U_asc = np.fromfile(U_file_asc, dtype='float32').reshape((length_asc, width_asc))

E_desc = np.fromfile(E_file_desc, dtype='float32').reshape((length_desc, width_desc))
N_desc = np.fromfile(N_file_desc, dtype='float32').reshape((length_desc, width_desc))
U_desc = np.fromfile(U_file_desc, dtype='float32').reshape((length_desc, width_desc))

# load the naf fault trace
fault_trace = np.loadtxt('data/naf_trace.xy')















# pre-allocate
vel_E = np.zeros((len(lat_regrid), len(lon_regrid)))
vel_U = np.zeros((len(lat_regrid), len(lon_regrid)))

# loop through every pixel
for ii in np.arange(0,len(lat_regrid)):
    for jj in np.arange(0,len(lon_regrid)):
        
        # create the design matrix
        G = np.array([[U_asc_regrid[ii,jj], E_asc_regrid[ii,jj]], [U_desc_regrid[ii,jj], E_desc_regrid[ii,jj]]])
        
        # get the two velocities for this pixel
        d = np.array([[vel_asc_regrid[ii,jj], vel_desc_regrid[ii,jj]]]).T
        
        # solve the linear system for the Up and East velocities
        m = np.linalg.solve(G, d)
        
        # save to arrays
        vel_U[ii,jj] = m[0]
        vel_E[ii,jj] = m[1]
        
