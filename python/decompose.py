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

#def decompose_xr(asc, desc, aschead, deschead, ascinc, descinc):
#    
U,E = decompose('Aok.tif', 'Dok.tif', aschead, deschead, ascinc, descinc)
aa = rioxarray.open_rasterio('Aok.tif')
aa.values[0]=U
export_xr2tif(aa,'U.tif', lonlat=False,dogdal=False)

import LiCSBAS_io_lib as io

def decompose(asctif, desctif, aschead, deschead, ascinc, descinc):
    vel_asc = io.read_geotiff(asctif)
    vel_desc = io.read_geotiff(desctif)
    vel_E = np.zeros(vel_desc.shape)
    vel_U = np.zeros(vel_desc.shape)
    #
    U_asc = np.cos(np.radians(ascinc))
    E_asc = -np.sin(np.radians(ascinc))*np.cos(np.radians(aschead))
    U_desc = np.cos(np.radians(descinc))
    E_desc = -np.sin(np.radians(descinc))*np.cos(np.radians(deschead))
    #
    for ii in np.arange(0,vel_E.shape[0]):
        for jj in np.arange(0,vel_E.shape[1]):
            # create the design matrix
            G = np.array([[U_asc, E_asc], [U_desc, E_desc]])
            # velocities
            d = np.array([[vel_asc[ii,jj], vel_desc[ii,jj]]]).T
            # solve the linear system for the Up and East velocities
            m = np.linalg.solve(G, d)
            # save to arrays
            vel_U[ii,jj] = m[0]
            vel_E[ii,jj] = m[1]
    return vel_U, vel_E


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
        
