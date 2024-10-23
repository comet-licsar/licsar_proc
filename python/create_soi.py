#!/usr/bin/env python3
'''
Created on 10/02/2024

@author: M. Nergizci, C. Magnard, M. Lazecky
University of Leeds
Gamma Remote Sensing
contact: mr.nergizci@gmail.com
'''
import numpy as np
import os
import py_gamma as pg
import sys
import subprocess
import shutil
import time
import matplotlib.pyplot as plt
from shapely.geometry import Polygon
from shapely.geometry import Point
from matplotlib.path import Path
from modules_sw_mn import * #functions saved here
import lics_unwrap as lu
import pickle
import LiCSAR_misc as misc



if len(sys.argv) < 2:
    print('Please provide pair information (assuming you run this in your BATCH_CACHE_DIR/frame)')
    print('e.g. python create_soi.py 20230129_20230210')
    sys.exit(1)

BLUE = '\033[94m'
ORANGE= '\033[38;5;208m'
ENDC = '\033[0m'  # ANSI code to end formatting

##Let's start!
start_time=time.time()

# print(BLUE + 'SBOVL is creating, please wait!!!' + ENDC)

##variables
tempdir = os.getcwd()
frame = os.path.basename(tempdir)
batchdir = os.environ['BATCH_CACHE_DIR'] #this is necesarry because the framebatch_gapfill run the code in $LiCS_temp folder which not correct location.
framedir = os.path.join(batchdir, frame)
pair=sys.argv[1]
#batchdir=os.environ['BATCH_CACHE_DIR']
prime, second = pair.split('_')
#framedir=os.path.join(batchdir, frame)
tr= int(frame[:3])
metafile = os.path.join(os.environ['LiCSAR_public'], str(tr), frame, 'metadata', 'metadata.txt')
master = misc.grep1line('master=',metafile).split('=')[1]
master_slcdir = os.path.join(framedir, 'SLC', master)
GEOC_folder=os.path.join(framedir, 'GEOC')
IFG_folder=os.path.join(framedir, 'IFG')
RSLC_folder=os.path.join(framedir, 'RSLC')
tab_folder=os.path.join(framedir, 'tab')

if (not os.path.exists(os.path.join(RSLC_folder, prime))) or (not os.path.exists(os.path.join(RSLC_folder, second))):
    print('ERROR: some of the input epoch data does not exist in your RSLC folder')
    sys.exit(1)


##output_dir
temp_file=os.path.join(framedir, 'temp_data')
if not os.path.exists(temp_file):
    os.makedirs(temp_file)


# File paths for stdout and stderr
stdout_log_path = os.path.join(temp_file, f"{pair}_soi_out.txt")
stderr_log_path = os.path.join(temp_file, f"{pair}_soi_err.txt")

# Open separate files for stdout and stderr
out_file = open(stdout_log_path, "w")
err_file = open(stderr_log_path, "w")

# Redirect stdout and stderr to their respective files
sys.stdout = out_file
sys.stderr = err_file

#tab_files
SLC1_tab_name= create_tab_file(prime, framedir, frame, type='RSLC')
RSLC2_tab_name= create_tab_file(second, framedir, frame, type='RSLC')
master_tab_name= create_tab_file(master, framedir, frame, type='SLC')
master_tab_name= create_tab_file(master, framedir, frame, type='RSLC')

#off
os.makedirs(os.path.join(IFG_folder, pair), exist_ok=True)
off_par = os.path.join(framedir, 'IFG', pair, pair + '.off')
# sim_unw = os.path.join(framedir, 'IFG', pair, pair+'.sim_unw') ##we don't need this because double differencing cancel out the topography correlated inf.
if not os.path.exists(off_par):  # Corrected os.file.exists to os.path.exists
    mpar = os.path.join(RSLC_folder, prime, prime + '.rslc.par')
    spar = os.path.join(RSLC_folder, second, second + '.rslc.par')
    
    # Command to create the offset
    exec_str = ['create_offset', mpar, spar, off_par, '1', '20', '4', '0']
    try:
        # Run the command and suppress the output
        subprocess.run(exec_str, check=True, stdout=subprocess.DEVNULL)
        # print(f"Command executed successfully: {' '.join(exec_str)}")
    except subprocess.CalledProcessError as e:
        print(f"An error occurred while executing the command: {e}")

SLC1_tab_mod1_name = os.path.join(tab_folder, prime+'_mod1.SLC_tab')
SLC1_tab_mod2_name = os.path.join(tab_folder, prime+'_mod2.SLC_tab')
RSLC2_tab_mod1_name = os.path.join(tab_folder, second+'_mod1.SLC_tab')
RSLC2_tab_mod2_name = os.path.join(tab_folder, second+'_mod2.SLC_tab')
SLC1_mod1_name = os.path.join(RSLC_folder,prime, prime+'_mod1.slc')
SLC1_mod2_name = os.path.join(RSLC_folder,prime, prime+'_mod2.slc')
RSLC2_mod1_name = os.path.join(RSLC_folder,second, second+'_mod1.slc')
RSLC2_mod2_name = os.path.join(RSLC_folder,second, second+'_mod2.slc')
mli1_mod1_name = os.path.join(RSLC_folder,prime, prime+'_mod1.mli')
mli1_mod2_name = os.path.join(RSLC_folder,prime, prime+'_mod2.mli')
rmli2_mod1_name = os.path.join(RSLC_folder,second, second+'_mod1.mli')
rmli2_mod2_name = os.path.join(RSLC_folder,second, second+'_mod2.mli')
diff_mod1_name = os.path.join(IFG_folder,pair, pair+'_mod1.diff')
diff_mod2_name = os.path.join(IFG_folder,pair, pair+'_mod2.diff')
diff_mod1_mask_name = os.path.join(IFG_folder,pair, pair+'_mod1_mask.diff')
diff_mod2_mask_name = os.path.join(IFG_folder,pair, pair+'_mod2_mask.diff')

#####variables for geocoding
lt_fine_suffix='lt_fine'
geo_dir= os.path.join(framedir, 'geo')
if os.path.exists(geo_dir) and os.path.isdir(geo_dir):
  for file in os.listdir(geo_dir):
    if file.endswith(lt_fine_suffix):
      lt_fine_file=os.path.join(geo_dir, file) 

  EQA_path=os.path.join(geo_dir, 'EQA.dem_par')
  EQA_par=pg.ParFile(EQA_path)
  widthgeo=EQA_par.get_value('width', dtype = int, index= 0)
  print(f'widthgeo value is {widthgeo}')

else:
  print(f'geo folder doesnt exists. Please check your {framedir}')
#####

print(BLUE + 'STEP-1: Double Difference Interferogram Generation!!!' + ENDC)

# number of looks
# here should be rechecked and make flexible for different multilooking. it should be improved!
rlks = 20
azlks = 4

# read SLC_tab files
SLC1_tab_temp = pg.read_tab(SLC1_tab_name)
RSLC2_tab_temp = pg.read_tab(RSLC2_tab_name)
master_tab_temp = pg.read_tab(master_tab_name)
SLC1_tab = np.array(framepath_tab(SLC1_tab_temp, framedir))
RSLC2_tab = np.array(framepath_tab(RSLC2_tab_temp, framedir))
master_tab=np.array(framepath_tab(master_tab_temp, framedir))

##create the path for SLC1_tab
SLC1_tab_long_name=os.path.join(framedir, 'tab', prime + '_tab_long')
master_tab_long_name=os.path.join(framedir, 'tab', master +'_tab_long')
# Open the file in write mode and write the formatted paths
if not os.path.exists(master_tab_long_name):
    with open(master_tab_long_name, 'w') as file:
        for line in master_tab:
            # Join the paths in the current line with a space
            combined_line = " ".join(line)
            # Write the combined line to the file
            file.write(combined_line + "\n")

if not os.path.exists(SLC1_tab_long_name):
    # Open the file in write mode and write the formatted paths
    with open(SLC1_tab_long_name, 'w') as file:
        for line in SLC1_tab:
            # Join the paths in the current line with a space
            combined_line = " ".join(line)
            # Write the combined line to the file
            file.write(combined_line + "\n")


nrows = SLC1_tab.shape[0]
ncols = SLC1_tab.shape[1]
# read image sizes
nr = []
naz = []

####Reading the un-multilooked and un-mosaic pixel number of range and azimuth from .slc.par file for each subswath.
for i in range(nrows):
  SLC_par = pg.ParFile(SLC1_tab[i][1])
  nr.append(SLC_par.get_value('range_samples', dtype = int, index = 0))
  naz.append(SLC_par.get_value('azimuth_lines', dtype = int, index = 0))

##prepare data with empty subswaths (create an array and store it as a binary file)
for i in range(nrows):
    # Create the file name using f-string
    file_name = f'empty.iw{1+i}.slc'
    full_path = os.path.join(temp_file, file_name)
    
    # Create the array
    if not os.path.exists(full_path):
        pg.create_array(full_path, nr[i], naz[i], 5, 0.0, 0.0) #output width nlines dtpes val val_im
      
SLC1_tab_mod1 = SLC1_tab.copy()
SLC1_tab_mod2 = SLC1_tab.copy()
RSLC2_tab_mod1 = RSLC2_tab.copy()
RSLC2_tab_mod2 = RSLC2_tab.copy()

##For both reference and secondary image, 2 SLC mosaics are generated: mod2 use subswaths 1 and 3, mod1 uses subswath 2. 
SLC1_tab_mod1[1][0] =os.path.join(temp_file, 'empty.iw2.slc')
SLC1_tab_mod2[0][0] =os.path.join(temp_file,'empty.iw1.slc')
SLC1_tab_mod2[2][0] =os.path.join(temp_file, 'empty.iw3.slc')
RSLC2_tab_mod1[1][0] =os.path.join(temp_file, 'empty.iw2.slc')
RSLC2_tab_mod2[0][0] =os.path.join(temp_file, 'empty.iw1.slc')
RSLC2_tab_mod2[2][0] =os.path.join(temp_file, 'empty.iw3.slc')

##save the new tabs to temp_file: mod1 is IW2, mod2 is IW1 and IW3
pg.write_tab(SLC1_tab_mod1, SLC1_tab_mod1_name)
pg.write_tab(SLC1_tab_mod2, SLC1_tab_mod2_name)
pg.write_tab(RSLC2_tab_mod1, RSLC2_tab_mod1_name)
pg.write_tab(RSLC2_tab_mod2, RSLC2_tab_mod2_name)


#replace burst_win range paramaters by ext_burst_win range paramaters:
for i in range(nrows):
    try:
        TOPS_par_master = pg.ParFile(master_tab[i][2])
        TOPS_par = pg.ParFile(SLC1_tab[i][2])
        nburst = TOPS_par.get_value('number_of_bursts', dtype=int, index=0)
        for b in range(nburst):
            ext_burst_win = TOPS_par_master.get_value(f'ext_burst_win_{b+1}')
            burst_win = TOPS_par.get_value(f'burst_win_{b+1}')

            TOPS_par.set_value(f'burst_win_{b+1}', ext_burst_win[0], index=0)
            TOPS_par.set_value(f'burst_win_{b+1}', ext_burst_win[1], index=1)
            TOPS_par.set_value(f'burst_win_{b+1}', ext_burst_win[4], index=4)
            TOPS_par.set_value(f'burst_win_{b+1}', ext_burst_win[5], index=5)

        TOPS_par.write_par(SLC1_tab[i][2])
        pg.update_par(SLC1_tab[i][2], SLC1_tab[i][2])

        ### Same process for RSLC
        TOPS_par = pg.ParFile(RSLC2_tab[i][2])
        nburst = TOPS_par.get_value('number_of_bursts', dtype=int, index=0)

        for b in range(nburst):
            ext_burst_win = TOPS_par_master.get_value(f'ext_burst_win_{b+1}')
            TOPS_par.set_value(f'burst_win_{b+1}', ext_burst_win[0], index=0)
            TOPS_par.set_value(f'burst_win_{b+1}', ext_burst_win[1], index=1)
            TOPS_par.set_value(f'burst_win_{b+1}', ext_burst_win[4], index=4)
            TOPS_par.set_value(f'burst_win_{b+1}', ext_burst_win[5], index=5)

        TOPS_par.write_par(RSLC2_tab[i][2])
        pg.update_par(RSLC2_tab[i][2], RSLC2_tab[i][2])

        ### Same process for master
        TOPS_par = pg.ParFile(master_tab[i][2])
        nburst = TOPS_par.get_value('number_of_bursts', dtype=int, index=0)

        for b in range(nburst):
            ext_burst_win = TOPS_par_master.get_value(f'ext_burst_win_{b+1}')
            TOPS_par.set_value(f'burst_win_{b+1}', ext_burst_win[0], index=0)
            TOPS_par.set_value(f'burst_win_{b+1}', ext_burst_win[1], index=1)
            TOPS_par.set_value(f'burst_win_{b+1}', ext_burst_win[4], index=4)
            TOPS_par.set_value(f'burst_win_{b+1}', ext_burst_win[5], index=5)

        TOPS_par.write_par(master_tab[i][2])
        pg.update_par(master_tab[i][2], master_tab[i][2])    

    except Exception as e:
        print(f"Error processing row {i}: {e}")


print(BLUE + 'mosaicking step!!' + ENDC)

# Check conditions and run pg.SLC_mosaic_ScanSAR if appropriate
pg.SLC_mosaic_ScanSAR(SLC1_tab_mod1_name, SLC1_mod1_name, SLC1_mod1_name + '.par', rlks, azlks, 0, master_tab_long_name)

pg.SLC_mosaic_ScanSAR(SLC1_tab_mod2_name, SLC1_mod2_name, SLC1_mod2_name + '.par', rlks, azlks, 0, master_tab_long_name)

pg.SLC_mosaic_ScanSAR(RSLC2_tab_mod1_name, RSLC2_mod1_name, RSLC2_mod1_name + '.par', rlks, azlks, 0, master_tab_long_name)

pg.SLC_mosaic_ScanSAR(RSLC2_tab_mod2_name, RSLC2_mod2_name, RSLC2_mod2_name + '.par', rlks, azlks, 0, master_tab_long_name)


print(BLUE + 'multilooking step!!' + ENDC)

##multilook
pg.multi_look(SLC1_mod1_name, SLC1_mod1_name + '.par', mli1_mod1_name, mli1_mod1_name + '.par', rlks, azlks)

pg.multi_look(SLC1_mod2_name, SLC1_mod2_name + '.par', mli1_mod2_name, mli1_mod2_name + '.par', rlks, azlks)

pg.multi_look(RSLC2_mod1_name, RSLC2_mod1_name + '.par', rmli2_mod1_name, rmli2_mod1_name + '.par', rlks, azlks)

pg.multi_look(RSLC2_mod2_name, RSLC2_mod2_name + '.par', rmli2_mod2_name, rmli2_mod2_name + '.par', rlks, azlks)

mli_mosaic_par = pg.ParFile(mli1_mod1_name + '.par')
mli_mosaic_nr = mli_mosaic_par.get_value('range_samples', dtype = int, index = 0)



# calculation of differential interferograms
# if not (os.path.exists(diff_mod1_name) and os.path.exists(diff_mod2_name)):
pg.SLC_intf(SLC1_mod1_name, RSLC2_mod1_name, SLC1_mod1_name + '.par', RSLC2_mod1_name + '.par', off_par, diff_mod1_name, rlks, azlks, 0, '-', 0, 0)
pg.SLC_intf(SLC1_mod2_name, RSLC2_mod2_name, SLC1_mod2_name + '.par', RSLC2_mod2_name + '.par', off_par, diff_mod2_name, rlks, azlks, 0, '-', 0, 0)

#sim_unw is not necesarry so I have used SLC_intf rather than SLC_diff_intf but it is example if I need in the future:
# pg.SLC_diff_intf(SLC1_mod2_name, RSLC2_mod2_name, SLC1_mod2_name + '.par', RSLC2_mod2_name + '.par', off_par, sim_unw, diff_mod2_name, rlks, azlks, 1, 0, 0.2, 1, 1)

##create bmp images to apply masking to extract overlaps
pg.rasmph_pwr(diff_mod1_name, mli1_mod1_name, mli_mosaic_nr)
pg.rasmph_pwr(diff_mod2_name, mli1_mod2_name, mli_mosaic_nr)

# mask data
pg.mask_data(diff_mod1_name, mli_mosaic_nr, diff_mod1_mask_name, diff_mod2_name + '.bmp', 1)
pg.mask_data(diff_mod2_name, mli_mosaic_nr, diff_mod2_mask_name, diff_mod1_name + '.bmp', 1)

##Redundant interval data, open if you need
# pg.rasmph_pwr(diff_mod1_mask_name, mli1_mod1_name, mli_mosaic_nr)
# pg.rasmph_pwr(diff_mod2_mask_name, mli1_mod2_name, mli_mosaic_nr)
# visualize differential phase of subswath overlap areas
#pg.dis2ras(diff_mod1_mask_name + '.bmp', diff_mod2_mask_name + '.bmp')
#pg.dis2mph(diff_mod1_mask_name, diff_mod2_mask_name, mli_mosaic_nr, mli_mosaic_nr)

####make the double differencing!!
print(BLUE + 'double difference interferogram are calculating!!!!' + ENDC)
print('printing type and shape of the inf data')

mli_par_path=os.path.join(mli1_mod1_name + '.par')

##the mli_par_path should be exist, othercase can't work properly!
try:
    if os.path.exists(mli_par_path):
        with open(mli_par_path, 'r') as mli_par:
            for line in mli_par:
                if line.startswith('range_samples'):
                    width = int(line.split(':')[-1].strip())
                elif line.startswith('azimuth_lines'):
                    az_line = int(line.split(':')[-1].strip())
    else:
        print(f'mli par does not exist. Please check the path: {mli_par_path}')
        sys.exit(1)

except Exception as e:
    print(f'An error occurred: {e}')
    sys.exit(1)
          
diff_mod1=np.fromfile(diff_mod1_mask_name, np.complex64).byteswap()
diff_mod2=np.fromfile(diff_mod2_mask_name, np.complex64).byteswap()
double_diff= diff_mod2*np.conj(diff_mod1)
diff_double_mask_temp=os.path.join(IFG_folder,pair, pair + '_soi_raw')
double_diff.byteswap().tofile(diff_double_mask_temp)

# #### Redundant interval data, open if you need
# diff_double_mask_pha_temp=os.path.join(IFG_folder,pair, pair + '_soi_raw_pha')
# ######
# pg.rasmph_pwr(diff_double_mask_temp, mli1_mod1_name, mli_mosaic_nr)
# ######
# ##extract phase values from interferogram.
# exec_str= ['cpx_to_real', diff_double_mask_temp, diff_double_mask_pha_temp, str(width), '4']
# try:
#   subprocess.run(exec_str, check=True, stdout=subprocess.DEVNULL)
#   # print(f"Command executed successfully: {' '.join(exec_str)}")
# except subprocess.CalledProcessError as e:
#   print(f"An error occurred while executing the command: {e}")
# phase_data = diff_double_mask_pha_temp
# geoc_file=os.path.join(GEOC_folder,pair, pair + '_soi_pha_raw.geo')
# exec_str=['geocode_back', phase_data, str(width), lt_fine_file, geoc_file, str(widthgeo), '0', '0', '0']
# try:
#   subprocess.run(exec_str, check=True, stdout=subprocess.DEVNULL)
#   # print(f"Command executed successfully: {' '.join(exec_str)}")
# except subprocess.CalledProcessError as e:
#   print(f"An error occurred while executing the command: {e}")
    
# geoc_tif=os.path.join(GEOC_folder,pair, pair + '_soi_pha_raw.geo.tif')
# exec_str=['data2geotiff', EQA_path, geoc_file,'2', geoc_tif, '0.0' ]
# try:
#   subprocess.run(exec_str, check=True, stdout=subprocess.DEVNULL)
#   # print(f"Command executed successfully: {' '.join(exec_str)}")
# except subprocess.CalledProcessError as e:
#   print(f"An error occurred while executing the command: {e}")
# ####

## The subswath interferogram calculated perfectly until here if you don't see any error in your screen. 
##So next step is to divide the subswath overlap into backward and forward as their precision is different.
print(BLUE + 'STEP-2: Polygons extraction...' + ENDC)
print(BLUE + 'bwr and fwr polygons are calculating!' + ENDC)

prime_RSLC=os.path.join(framedir, 'RSLC', prime)
SLC_par=pg.ParFile(os.path.join(prime_RSLC, prime + '.rslc.par'))

nsubswaths = 3
TOPS_par = []
for i in range(1,4):
    TOPS_file=os.path.join(prime_RSLC, prime + f'.IW{i}.rslc.TOPS_par')
    if os.path.exists(TOPS_file):
        TOPS_par.append(pg.ParFile(TOPS_file))
            
# # read SLC parameters
r0 = SLC_par.get_value('near_range_slc', index = 0, dtype = float)        # near range distance
rps = SLC_par.get_value('range_pixel_spacing', index = 0, dtype = float)  # range pixel spacing
t0 = SLC_par.get_value('start_time', index = 0, dtype = float)            # time of first azimuth line
tazi = SLC_par.get_value('azimuth_line_time', index = 0, dtype = float)   # time between each azimuth line

# get number of bursts per subswath
# max_nbursts = 0
# nbursts = []
# for sw in range(nsubswaths):
#   nbursts.append(TOPS_par[sw].get_value('number_of_bursts', index = 0, dtype = int))
#   if nbursts[sw] > max_nbursts:
#     max_nbursts = nbursts[sw]

max_nbursts = 0
min_nbursts = float('inf')  # Initialize min_nbursts to infinity
nbursts = []

# Find minimum burst number, then fix for all subswaths.
for sw in range(nsubswaths):
    nb = TOPS_par[sw].get_value('number_of_bursts', index=0, dtype=int)  # Mocked function call
    print(f'number of burst in reality = {nb}')
    if nb > max_nbursts:
        max_nbursts = nb
    if nb < min_nbursts:
        min_nbursts = nb

# Update nbursts to minimum burst size across all subswaths
print(f'all burst number set to be {max_nbursts} in each subswath!')
nbursts = [max_nbursts] * nsubswaths

# initialize first and last range and azimuth pixel of each burst (SLC mosaic)
rpix0 = np.zeros((nsubswaths, max_nbursts))
rpix2 = np.zeros((nsubswaths, max_nbursts))
azpix0 = np.zeros((nsubswaths, max_nbursts))
azpix2 = np.zeros((nsubswaths, max_nbursts))

for sw in range(nsubswaths):
  for b in range(nbursts[sw]):
    # read burst window (modified burst window as in previous e-mail)
    ext_burst_win = TOPS_par[sw].get_value('ext_burst_win_%d' %(b+1))
    burst_win = TOPS_par[sw].get_value('burst_win_%d' %(b+1))
    # calculate pixel coordinates of bursts in mosaicked image
    if burst_win is not None:
        rpix0[sw, b] = round((float(burst_win[0]) - r0) / rps)  ## When you change the ext_burst_win overlap in azimuth direction become overtook.
        rpix2[sw, b] = round((float(burst_win[1]) - r0) / rps)
        azpix0[sw, b] = round((float(burst_win[2]) - t0) / tazi)
        azpix2[sw, b] = round((float(burst_win[3]) - t0) / tazi)
    else:
        rpix0[sw, b] = np.nan
        rpix2[sw, b] = np.nan
        azpix0[sw, b] = np.nan
        azpix2[sw, b] = np.nan

# First and last range and azimuth pixel of each burst (MLI mosaic / interferogram geometry)
rpix_ml0 = np.where(np.isnan(rpix0), np.nan, np.round(rpix0 / rlks).astype(int))
rpix_ml2 = np.where(np.isnan(rpix2), np.nan, np.round((rpix2 + 1) / rlks - 1).astype(int))
azpix_ml0 = np.where(np.isnan(azpix0), np.nan, np.round(azpix0 / azlks).astype(int))
azpix_ml2 = np.where(np.isnan(azpix2), np.nan, np.round((azpix2 + 1) / azlks - 1).astype(int))

# calculate intersection of bursts (subswath intersection)
p_inter = []
bursts = []
for sw in range(nsubswaths - 1):
    p_inter_sw = []
    for b0 in range(nbursts[sw]):
        if np.isnan(rpix_ml0[sw, b0]) or np.isnan(azpix_ml0[sw, b0]):
            p_inter_sw.append([None] * nbursts[sw+1])
            continue
        p_inter_b0 = []
        rg_az1 = [
            [rpix_ml0[sw, b0], azpix_ml0[sw, b0]],
            [rpix_ml2[sw, b0], azpix_ml0[sw, b0]],
            [rpix_ml2[sw, b0], azpix_ml2[sw, b0]],
            [rpix_ml0[sw, b0], azpix_ml2[sw, b0]],
            [rpix_ml0[sw, b0], azpix_ml0[sw, b0]],  # Close the polygon
        ]
        p0 = Polygon(rg_az1)
        bursts.append(p0)
        for b1 in range(nbursts[sw+1]):
            if np.isnan(rpix_ml0[sw+1, b1]) or np.isnan(azpix_ml0[sw+1, b1]):
                p_inter_b0.append(None)
                continue
            rg_az2 = [
                [rpix_ml0[sw+1, b1], azpix_ml0[sw+1, b1]],
                [rpix_ml2[sw+1, b1], azpix_ml0[sw+1, b1]],
                [rpix_ml2[sw+1, b1], azpix_ml2[sw+1, b1]],
                [rpix_ml0[sw+1, b1], azpix_ml2[sw+1, b1]],
                [rpix_ml0[sw+1, b1], azpix_ml0[sw+1, b1]],  # Close the polygon
            ]
            p1 = Polygon(rg_az2)
            p_inter_b0.append(p0.intersection(p1))
        p_inter_sw.append(p_inter_b0)
    p_inter.append(p_inter_sw)

bwr = []
fwr = []

counter = 0  # Initialize a counter to alternate between fwr and bwr

# Define the backward and forward subswath overlap polygons
for sw in range(len(p_inter)):
    # Reset counter when switching to a new sub swath
    if sw == 1:
        counter = 0
        
    # Iterate over bursts in the current sub swath
    for b0 in range(len(p_inter[sw])):
        for b1 in range(len(p_inter[sw][b0])):
            if p_inter[sw][b0][b1] is not None:
                exterior = p_inter[sw][b0][b1].exterior
                if exterior:
                    x, y = exterior.xy
                    if len(x) > 0 and len(y) > 0:
                        if (azpix_ml0[1][0] - azpix_ml0[0][0] > 0 and sw == 0):
                            if counter % 2 == 0:
                                fwr.append(p_inter[sw][b0][b1])
                            else:
                                bwr.append(p_inter[sw][b0][b1])
                            counter += 1
                        elif (azpix_ml0[1][0] - azpix_ml0[0][0] < 0 and sw == 0):
                            if counter % 2 == 0:
                                bwr.append(p_inter[sw][b0][b1])
                            else:
                                fwr.append(p_inter[sw][b0][b1])
                            counter += 1                            
                        elif (azpix_ml0[2][0] - azpix_ml0[1][0] > 0 and sw == 1):
                            if counter % 2 == 0:
                                fwr.append(p_inter[sw][b0][b1])
                            else:
                                bwr.append(p_inter[sw][b0][b1])
                            counter += 1
                        elif (azpix_ml0[2][0] - azpix_ml0[1][0] < 0 and sw == 1):
                            if counter % 2 == 0:
                                bwr.append(p_inter[sw][b0][b1])
                            else:
                                fwr.append(p_inter[sw][b0][b1])
                            counter += 1

##saving bwr and fwr as pickle##
bwrs=os.path.join(IFG_folder,pair,pair +'_bwr.pkl')
fwrs=os.path.join(IFG_folder,pair,pair + '_fwr.pkl')
with open(bwrs, 'wb') as f:
    pickle.dump(bwr, f)

# Saving fwr to a file named 'fwr.pkl'
with open(fwrs, 'wb') as f:
    pickle.dump(fwr, f)

print(BLUE + 'STEP-3:Scaling factor and filtering!' + ENDC)
#######calculating scaling factors
path_to_slcdir = os.path.join(framedir, 'RSLC', prime)
sf_array=get_sf_array(path_to_slcdir, f0=5405000500, burst_interval=2.758277)

###reopen interferogram
dd_cpx=np.fromfile(diff_double_mask_temp, np.complex64).byteswap().reshape(az_line, width)

print(BLUE + 'adf_filtering!!!!' + ENDC)
def process_diff_data(dd_cpx, polygons_filtered,sf_array, pair, suffix, width, az_line):
    poly_mask=generate_mask_for_polygons_optimized(dd_cpx.shape, polygons_filtered)
    poly_raw=np.where(poly_mask, dd_cpx, np.nan)
    # Replace NaN in real and imaginary parts with 0
    real_part_nonan = np.nan_to_num(poly_raw.real)
    imaginary_part_nonan = np.nan_to_num(poly_raw.imag)
    poly_raw_nonan = real_part_nonan + 1j * imaginary_part_nonan
    poly_raw_nonan_path = os.path.join(IFG_folder,pair, f"{pair}_{suffix}_soi_nonan")
    poly_raw_nonan.astype(np.complex64).byteswap().tofile(poly_raw_nonan_path)

    # ADF filtering
    poly_adf_path = os.path.join(IFG_folder,pair, f"{pair}_{suffix}_soi_adf")
    poly_pha_path = os.path.join(IFG_folder,pair, f"{pair}_{suffix}_soi_adf_pha")
    poly_coh_path = os.path.join(IFG_folder,pair, f"{pair}_{suffix}_soi_adf_coh")

    # Execute ADF
    exec_str = ['adf', poly_raw_nonan_path, poly_adf_path, poly_coh_path, str(width), '1', '-', '-', '-', '-', '-', '-']
    try:
        subprocess.run(exec_str, check=True, stdout=subprocess.DEVNULL)
        # print(f"ADF executed successfully for {suffix}")
    except subprocess.CalledProcessError as e:
        print(f"An error occurred while executing ADF for {suffix}: {e}")

    # Convert complex to real
    exec_str = ['cpx_to_real', poly_adf_path, poly_pha_path, str(width), '4']
    try:
        subprocess.run(exec_str, check=True, stdout=subprocess.DEVNULL)
        # print(f"Conversion to real executed successfully for {suffix}")
    except subprocess.CalledProcessError as e:
        print(f"An error occurred during conversion to real for {suffix}: {e}")

    # Load the phase data
    poly_pha = np.fromfile(poly_pha_path, dtype=np.float32).byteswap().reshape(az_line, width)
    print(f"Phase data loaded for {suffix}")

    ###scaling_factor:
    poly_pha_scaled=poly_pha*sf_array
    poly_pha_scaled_path = os.path.join(IFG_folder,pair, f"{pair}_{suffix}_soi_adf_scaled")
    poly_pha_scaled.astype(np.float32).byteswap().tofile(poly_pha_scaled_path)
    
    return

# Loop over both fwr and bwr data
for suffix, polygons_filtered in [('fwr', fwr), ('bwr', bwr)]:
    # Process each and store the resulting phase data
    process_diff_data(dd_cpx, polygons_filtered,sf_array, pair, suffix, width, az_line)

print('HouseKeeping')

#######################################
#####producing merged adf and cc values
bwr_ddif_coh_path=os.path.join(IFG_folder,pair, f"{pair}_bwr_soi_adf_coh")
fwr_ddif_coh_path=os.path.join(IFG_folder,pair, f"{pair}_fwr_soi_adf_coh")
bwr_ddif_coh=np.fromfile(bwr_ddif_coh_path, dtype=np.float32).byteswap().reshape(az_line, width)
fwr_ddif_coh=np.fromfile(fwr_ddif_coh_path, dtype=np.float32).byteswap().reshape(az_line, width)
##replace zero2nan
bwr_ddif_coh[bwr_ddif_coh==0]=np.nan
fwr_ddif_coh[fwr_ddif_coh==0]=np.nan

#merging coh..
ddif_coh_merg=np.copy(fwr_ddif_coh)
mask_bwr=~np.isnan(bwr_ddif_coh)
ddif_coh_merg[mask_bwr]=bwr_ddif_coh[mask_bwr]
ddif_coh_merg_path=os.path.join(IFG_folder,pair, f'{pair}_soi_adf_coh')
ddif_coh_merg.astype(np.float32).byteswap().tofile(ddif_coh_merg_path)

##same for merged adf values:
bwr_ddif_adf_path=os.path.join(IFG_folder,pair, f'{pair}_bwr_soi_adf_scaled')
fwr_ddif_adf_path=os.path.join(IFG_folder,pair, f'{pair}_fwr_soi_adf_scaled')
bwr_ddif_adf=np.fromfile(bwr_ddif_adf_path, dtype=np.float32).byteswap().reshape(az_line, width)
fwr_ddif_adf=np.fromfile(fwr_ddif_adf_path, dtype=np.float32).byteswap().reshape(az_line, width)

##replace zero2nan
bwr_ddif_adf[bwr_ddif_adf==0]=np.nan
fwr_ddif_adf[fwr_ddif_adf==0]=np.nan

#merging adf..
ddif_adf_merg=np.copy(fwr_ddif_adf)
mask_bwr=~np.isnan(bwr_ddif_adf)
ddif_adf_merg[mask_bwr]=bwr_ddif_adf[mask_bwr]
ddif_adf_merg=ddif_adf_merg.astype(np.float32)
ddiff_adf_merg_path=os.path.join(IFG_folder,pair, f'{pair}_soi_adf_scaled')
ddif_adf_merg.byteswap().tofile(ddiff_adf_merg_path)

##extra for only adf_nonscaled
##same for merged adf values:
bwr_ddif_adf_nonscale_path=os.path.join(IFG_folder,pair, f'{pair}_bwr_soi_adf_pha')
fwr_ddif_adf_nonscale_path=os.path.join(IFG_folder,pair, f'{pair}_fwr_soi_adf_pha')
bwr_ddif_adf_nonscale=np.fromfile(bwr_ddif_adf_nonscale_path, dtype=np.float32).byteswap().reshape(az_line, width)
fwr_ddif_adf_nonscale=np.fromfile(fwr_ddif_adf_nonscale_path, dtype=np.float32).byteswap().reshape(az_line, width)

##replace zero2nan
bwr_ddif_adf_nonscale[bwr_ddif_adf_nonscale==0]=np.nan
fwr_ddif_adf_nonscale[fwr_ddif_adf_nonscale==0]=np.nan

#merging adf..
ddif_adf_nonscale_merg=np.copy(fwr_ddif_adf_nonscale)
mask_bwr=~np.isnan(bwr_ddif_adf_nonscale)
ddif_adf_nonscale_merg[mask_bwr]=bwr_ddif_adf_nonscale[mask_bwr]
ddif_adf_nonscale_merg=ddif_adf_nonscale_merg.astype(np.float32)
ddiff_adf_nonscale_merg_path=os.path.join(IFG_folder,pair, f'{pair}_soi_adf_unscaled')
ddif_adf_nonscale_merg.byteswap().tofile(ddiff_adf_nonscale_merg_path)

for prefix in ['coh', 'scaled']:  #unscaled
    phase_data = os.path.join(IFG_folder,pair, f'{pair}_soi_adf_{prefix}')
    geoc_file=os.path.join(GEOC_folder,pair, f'{pair}_soi_adf_{prefix}.geo')
    exec_str=['geocode_back', phase_data, str(width), lt_fine_file, geoc_file, str(widthgeo), '0', '0', '0']
    try:
      subprocess.run(exec_str, check=True, stdout=subprocess.DEVNULL)
      # print(f"Command executed successfully: {' '.join(exec_str)}")
    except subprocess.CalledProcessError as e:
      print(f"An error occurred while executing the command: {e}")
        
    geoc_tif=os.path.join(GEOC_folder,pair, f'{pair}_soi_adf_{prefix}.geo.tif')
    exec_str=['data2geotiff', EQA_path, geoc_file,'2', geoc_tif, '0.0']
    try:
      subprocess.run(exec_str, check=True, stdout=subprocess.DEVNULL)
      # print(f"Command executed successfully: {' '.join(exec_str)}")
    except subprocess.CalledProcessError as e:
      print(f"An error occurred while executing the command: {e}")

# After script finishes, restore stdout and stderr
sys.stdout = sys.__stdout__
sys.stderr = sys.__stderr__

# Close the log files
out_file.close()
err_file.close()

end_time=time.time()
execution_time = end_time - start_time
print(f"SBOI generated in {round(execution_time)} seconds, You're welcome..")
