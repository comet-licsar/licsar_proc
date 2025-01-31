#!/usr/bin/env python3
# these are functions for decomposition into E,U vectors from 2 or more tracks

import subprocess as subp
import numpy as np
from scipy import interpolate
from lics_unwrap import *
import dask.array as da



'''
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
'''

def decompose_framencs(framencs, extract_cum = False, medianfix = False, annual = False, annual_buffer_months = 0):
    """ will decompose frame licsbas results
    the basenames in framencs should contain frame id, followed by '.', e.g.:
    framencs = ['062D_07629_131313.nc', '172A_07686_131012.nc']
    
    Args:
        framencs (list):  licsbas nc result files, named by their frame id
        extract_cum (bool): if True, will use the first frame and convert to pseudo vertical
        annual (bool):   if True, will estimate and decompose annual velocities
        annual_buffer_months (int): adds extra months for annual velocities
    Returns:
        xr.Dataset with U, E, [cum_vert] arrays
    """
    framesetvel = []
    frameset = []
    firstrun = True
    # getting years in all ncs:
    yearsall = None
    for nc in framencs:
        frame = os.path.basename(nc).split('.')[0]
        framenc = xr.open_dataset(nc)
        years = framenc.time.dt.year.values
        years = list(set(years))
        # print(years)
        if not yearsall:
            yearsall = years
        else:
            for y in yearsall.copy():
                if y not in years:
                    yearsall.remove(y)
    if len(yearsall)==0:
        print('no overlapping year, cancelling annual decomposition')
        return False
    for nc in framencs:
        frame = os.path.basename(nc).split('.')[0]
        print('extracting frame '+frame)
        inc, heading = get_frame_inc_heading(frame)
        framenc = xr.open_dataset(nc)
        framevel = framenc['vel']
        if medianfix:
            framevel = framevel - framevel.median()
        if firstrun:
            template = framevel.copy()
            firstrun = False 
            if extract_cum:
                cum_vert = framenc['cum']
                inc = inc.interp_like(framevel)
                cum_vert = cum_vert/np.cos(np.radians(inc))
        else:
            framevel = framevel.interp_like(template)
            if annual:
                framenc = framenc.interp_like(template)
        inc = inc.interp_like(framevel)
        heading = heading.interp_like(framevel)
        framesetvel.append((framevel.values, heading.values, inc.values))
        if annual:
            # doing the annuals!
            nc1 = calculate_annual_vels(framenc, yearsall, annual_buffer_months)
            frameset.append((nc1['vel_annual'], heading.values, inc.values))
    dec = xr.Dataset()
    U = template.copy()
    E = template.copy()
    U.values, E.values = decompose_np_multi(framesetvel)
    dec['U'] = U
    dec['E'] = E
    if annual:
        # if annual, then frameset is from nc.vel_annual, heading.values, inc.values
        years = None
        for framedata in frameset:
            if not years:
                years = list(framedata[0].year.values)
            else:
                yearst = list(framedata[0].year.values)
                for year in years:
                    if year not in yearst:
                        years.remove(year)
        # ok, now then decompose the annuals per year
        vUxr = frameset[0][0].sel(year = years).copy().rename('vU')*0
        vExr = frameset[0][0].sel(year = years).copy().rename('vE')*0
        for year in years:
            annualset = []
            for framedata in frameset:
                vel = framedata[0].sel(year = year).values
                annualset.append((vel, framedata[1], framedata[2]))
            print('decomposing year '+str(year))
            try:
                vU, vE = decompose_np_multi(annualset) #, beta = 0)
                vUxr.loc[year,:,:] = vU
                vExr.loc[year,:,:] = vE
            except:
                print('error decomposing, setting nans')
        dec['vU'] = vUxr
        dec['vE'] = vExr
    if extract_cum:
        dec['cum'] = cum_vert
    return dec



import LiCSBAS_inv_lib as inv_lib

# Copilot-generated function for custom resample:
import pandas as pd
# cube is either xr.Dataset or xr.Dataarray
def custom_annual_resample(cube, buffermonths=6):
    start_date = pd.to_datetime(cube.time.values[0])
    end_date = pd.to_datetime(cube.time.values[-1])
    #
    # Generate a new time range that extends 6 months before and after
    new_time_range = pd.date_range(
        start=start_date - pd.DateOffset(months=buffermonths),
        end=end_date + pd.DateOffset(months=buffermonths),
        freq='M'
    )
    #
    # Create an empty list to hold the new resampled cubes
    resampled_cubes = []
    #
    for date in pd.date_range(start=start_date, end=end_date, freq='AS'):
        # Define the range for the current resample
        start_period = date - pd.DateOffset(months=buffermonths)
        end_period = date + pd.DateOffset(months=12+buffermonths) - pd.DateOffset(days=1)
        #
        # Select the data within this range
        selected_data = cube.sel(time=slice(start_period, end_period))
        #
        # Append the selected data to the list
        resampled_cubes.append(selected_data)
    #
    # Combine all resampled data into one DataArray or Dataset
    resampled_cube = xr.concat(resampled_cubes, dim='time')
    #
    return resampled_cube


def calculate_annual_vels(cube, commonyears = None, buffermonths = 0):
    """Will calculate annual velocities from LiCSBAS results
    
    Args:
        cube (xr.Dataset): loaded netcdf file, extracted using e.g. LiCSBAS_out2nc.py
        commonyears (list): list of years to decompose
        buffermonths (int): extend selection of annual data by a number of +-buffermonths months (experimental)
    Returns:
        xr.Dataset with new dataarray: vel_annual
    """
    if commonyears:
        cube = cube.sel(time=np.isin(cube.time.dt.year.values, commonyears))
    if buffermonths > 0:
        print('Warning, using Copilot-generated trick to add more months around year of interest...')
        annualset = custom_annual_resample(cube, buffermonths)
    else:
        annualset = cube.cum.resample(time='AS')
    firstrun = True
    for yeardt, yearcum in annualset:
        year = str(yeardt).split('-')[0]
        print('processing year '+year)
        dt_cum = (np.array([tdate.toordinal() for tdate in yearcum.time.dt.date.values]) - pd.Timestamp(yeardt).date().toordinal())/365.25 # in fraction of a year
        # see LiCSBAS_cum2vel.py:
        cum_tmp = yearcum.values
        n_im, length, width = cum_tmp.shape
        bool_allnan = np.all(np.isnan(cum_tmp), axis=0)
        vconst = np.zeros((length, width), dtype=np.float32)*np.nan
        vel = np.zeros((length, width), dtype=np.float32)*np.nan
        #
        cum_tmp = cum_tmp.reshape(n_im, length*width)[:, ~bool_allnan.ravel()].transpose()
        vel[~bool_allnan], vconst[~bool_allnan] = inv_lib.calc_vel(cum_tmp, dt_cum)
        #vel[~bool_allnan], vconst[~bool_allnan] = calc_vel(cum_tmp, dt_cum)
        vel_annual = cube.vel.copy()
        vel_annual.values = vel
        vel_annual = vel_annual.assign_coords({'year':int(year)}).expand_dims('year').rename('vel_annual')
        if firstrun:
            vel_annual_cube = vel_annual.copy()
            firstrun = False
            #dv = vel
        else:
            vel_annual_cube = xr.concat([vel_annual_cube,vel_annual], dim='year')
            #dv = np.vstack(dv,vel)
    cube['vel_annual'] = vel_annual_cube.copy()
    return cube




def get_frame_inc_heading(frame):
    """Extracts inc and heading from E, U for given frame
    """
    geoframedir = os.path.join(os.environ['LiCSAR_public'], str(int(frame[:3])), frame)
    # look angle (inc) / heading - probably ok, but needs check:
    e=os.path.join(geoframedir,'metadata',frame+'.geo.E.tif')
    #n=os.path.join(geoframedir,'metadata',frame+'.geo.N.tif') #no need for N
    u=os.path.join(geoframedir,'metadata',frame+'.geo.U.tif')
    return extract_inc_heading(e, u)


def extract_inc_heading(efile, ufile):
    e = load_tif2xr(efile)
    e = e.where(e != 0)
    #n = load_tif2xr(n, cliparea_geo=cliparea)
    u = load_tif2xr(ufile)
    u = u.where(u != 0)
    #
    theta=np.arcsin(u)
    phi=np.arccos(e/np.cos(theta))
    heading = np.rad2deg(phi)-180
    inc = 90-np.rad2deg(theta)   #correct
    #inc.values.tofile(outinc)
    return inc, heading


def decompose_dask(cube, blocklen=5, num_workers=5):
    """Simple parallel decomposition of dec. datacube (must have asc,desc,asc_inc,desc_inc, asc_heading, desc_heading arrays)
    """
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


def decompose_xr(asc, desc, heading_asc, heading_desc, inc_asc, inc_desc, beta=0):
    """Perform simple decomposition for two frames in asc and desc.
    inputs are xr.dataarrays - this will also check/interpolate them to fit
    """
    cube=xr.Dataset()
    cube['asc'] = asc
    cube['desc'] = desc.interp_like(asc, method='linear'); desc=None
    cube['U']=cube.asc.copy()
    cube['E']=cube.asc.copy()
    if not np.isscalar(heading_asc):
        cube['asc_heading'] = heading_asc.interp_like(asc, method='linear'); heading_asc=cube.asc_heading.values
        cube['desc_heading'] = heading_desc.interp_like(asc, method='linear'); heading_desc=cube.desc_heading.values
    if not np.isscalar(inc_asc):
        cube['asc_inc'] = inc_asc.interp_like(asc, method='linear'); inc_asc=cube.asc_inc.values
        cube['desc_inc'] = inc_desc.interp_like(asc, method='linear'); inc_desc=cube.desc_inc.values
    cube['U'].values, cube['E'].values = decompose_np(cube.asc.values, cube.desc.values, heading_asc, heading_desc, inc_asc , inc_desc)
    return cube[['U', 'E']]


# 2022-10-18 - this should be pretty good one (next only use weights or something)
def decompose_np(vel_asc, vel_desc, aschead, deschead, ascinc, descinc, beta=0, do_velUN = True):
    """Decomposes values from ascending and descending np (or xr) arrays, using heading and inc. angle
    (these might be arrays of same size of just float values)
    
    Args:
        beta (float): angle of expected horizontal motion direction, clockwise from the E, in degrees
        do_velUN (boolean): if yes, extract v_{UN} instead of v_U, following Qi et al, 2022: doi=10.1029/2F2022JB024176
    """
    vel_E = np.zeros(vel_desc.shape)
    vel_U = np.zeros(vel_desc.shape)
    #
    if do_velUN:
        U_asc = np.sqrt(1 - (np.sin(np.radians(ascinc))**2) * (np.cos(np.radians(aschead))**2))
        U_desc = np.sqrt(1 - (np.sin(np.radians(descinc)) ** 2) * (np.cos(np.radians(deschead)) ** 2))
    else:
        U_asc = np.cos(np.radians(ascinc))
        U_desc = np.cos(np.radians(descinc))
    #
    E_asc = -np.sin(np.radians(ascinc))*np.cos(np.radians(aschead+beta))
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
this is to load 3 datasets and decompose them:
dirpath='/gws/nopw/j04/nceo_geohazards_vol1/public/shared/temp/earmla'
#for frame in []
nc1 = os.path.join(dirpath, '051D_03973_131313.nc')
nc1=xr.open_dataset(nc1)
vel1 = nc1.vel.values
heading1 = -169.87
inc1 = 43.64

nc2 = os.path.join(dirpath, '124D_04017_131313.nc')
nc2=xr.open_dataset(nc2)
vel2 = nc2.vel.interp_like(nc1.vel).values
heading2 = -169.88
inc2 = 34.98

nc3 = os.path.join(dirpath, '175A_03997_131313.nc')
nc3=xr.open_dataset(nc3)
vel3 = nc3.vel.interp_like(nc1.vel).values
heading3 = -10.16
inc3 = 38.42

years = np.array([2015, 2016, 2017, 2018, 2019, 2020, 2021, 2022])
vUxr = nc1.vel_annual.sel(year = years).copy().rename('vU')
vExr = nc1.vel_annual.sel(year = years).copy().rename('vE')
for year in years:
    vel1 = nc1.vel_annual.sel(year = year).values
    vel2 = nc2.vel_annual.sel(year = year).values
    vel3 = nc3.vel_annual.sel(year = year).values
    input_data = [(vel1, heading1, inc1), (vel2, heading2, inc2), (vel3, heading3, inc3)]
    print('decomposing year '+str(year))
    vU, vE = decompose_np_multi(input_data, beta = 0)
    vUxr.loc[year,:,:] = vU
    vExr.loc[year,:,:] = vE


decomposedxr = xr.Dataset()
decomposedxr['vU'] = vUxr
decomposedxr['vE'] = vExr
decomposedxr.to_netcdf('decomposed_s1.nc')
'''


def decompose_np_multi(input_data, beta = 0):
    """Decompose 2 or more frames
    
    Args:
        input data (list of tuples) e.g. input_data = [(vel1, heading1, inc1), (vel2, heading2, inc2), (vel3, heading3, inc3)]
    
    Note: velX is np.array and headingX/incX is in degrees, either a number or np.array
    """
    #
    template = input_data[0][0]
    vel_E = np.zeros(template.shape)
    vel_U = np.zeros(template.shape)
    #
    Us=list()
    Es=list()
    vels = []
    for frame in input_data:
        vel = frame[0]
        heading = frame[1]
        incangle = frame[2]
        #Us = np.append(Us, np.cos(np.radians(incangle)))
        Us.append(np.cos(np.radians(incangle)))
        #Es = np.append(Es, -np.sin(np.radians(incangle))*np.cos(np.radians(heading+beta)))
        Es.append(-np.sin(np.radians(incangle))*np.cos(np.radians(heading+beta)))
        vels.append(vel)
        # run for each pixel
    numframes = len(vels)
    Us = np.array(Us)
    Es = np.array(Es)
    for ii in np.arange(0,vel_E.shape[0]):
        for jj in np.arange(0,vel_E.shape[1]):
            # prepare template for d = G m
            d = np.array(())
            for i in range(numframes):
                d = np.append(d, np.array([vels[i][ii,jj]]))
            d = np.array([d]).T
            if np.isnan(d).any():
                # if at least one is nan, skip it:  # can improve it but 'all' is not an option
                #if np.isnan(np.max(d)):
                vel_U[ii,jj] = np.nan
                vel_E[ii,jj] = np.nan
            else:
                # create the design matrix
                if np.isscalar(Us[0]):  # in case of only values (i.e. one inc and heading per each frame)
                    G = np.vstack([Us, Es]).T              
                else:  # in case this is array  # not tested!
                    G = np.vstack([Us[:,ii,jj], Es[:,ii,jj]]).T
                # solve the linear system for the Up and East velocities
                #m = np.linalg.solve(G, d)
                try:
                    m = np.linalg.lstsq(G, d)[0]
                    # save to arrays
                    vel_U[ii,jj] = m[0]
                    vel_E[ii,jj] = m[1]
                except:
                    vel_U[ii,jj] = np.nan
                    vel_E[ii,jj] = np.nan
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


''' (old) usage example:
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

'''




'''
orig AW approach:
import matplotlib.pyplot as plt


# these packages are only needed for the final multivariate plot
import seaborn as sns
import pandas as pd
import interseis_lib as lib



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
        
'''
