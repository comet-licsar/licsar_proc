#!/usr/bin/env python3
'''
Last Revision August 2024

@author: M. Nergizci
University of Leeds,
contact: mr.nergizci@gmail.com
'''
import numpy as np
import os
import sys
import subprocess
import shutil
from modules_sw_mn import * #functions saved here


if len(sys.argv) < 2:
    print('Please provide pair information: i.e python create_sbovl_ifg.py 20230129_20230210')
    sys.exit(1)

# User feedback with ANSI colors
BLUE = '\033[94m'
ORANGE = '\033[38;5;208m'
ENDC = '\033[0m'  # ANSI code to end formatting

##variables
framedir = os.getcwd()
frame = os.path.basename(os.getcwd())
#TODO I can check the frame name is okay for format of LiCSAR_frame.
pair = sys.argv[1]
#batchdir = os.environ['BATCH_CACHE_DIR']
prime, second = pair.split('_')
#framedir = os.path.join(batchdir, frame)
GEOC_folder = os.path.join(framedir, 'GEOC')

# Define the paths
boi_path = os.path.join(GEOC_folder, pair, pair + '.geo.bovldiff.adf.mm.tif')
soi_path = os.path.join(GEOC_folder, pair, pair + '_soi_adf_scaled.geo.tif')
output_path = os.path.join(GEOC_folder, pair, pair + '.geo.sbovldiff.adf.mm.tif')

# Coherence paths
boi_coh_path = os.path.join(GEOC_folder, pair, pair + '.geo.bovldiff.adf.cc.tif')
soi_coh_path = os.path.join(GEOC_folder, pair, pair + '_soi_adf_coh.geo.tif')
output_coh_path = os.path.join(GEOC_folder, pair, pair + '.geo.sbovldiff.adf.cc.tif')

# Open input GeoTIFF files
if os.path.exists(boi_path):
    boi_mm = open_geotiff(boi_path, fill_value=np.nan)
else:
    print(f"{boi_path} doesn't exist ,please run create_bovl_ifg.sh first!")
    sys.exit(1)
    
if os.path.exists(soi_path):
    soi = open_geotiff(soi_path, fill_value=np.nan)
else:
    print(f"{soi_path} doesn't exist ,please run create_soi.py first!")
    sys.exit(1)

# Open coherences
if os.path.exists(boi_coh_path):
    boi_coh = open_geotiff(boi_coh_path, fill_value=np.nan)
else:
    print(f"{boi_coh_path} doesn't exist ,please run create_bovl_ifg.sh first!")
    sys.exit(1)
if os.path.exists(soi_coh_path):
    soi_coh = open_geotiff(soi_coh_path, fill_value=np.nan)
else:
    print(f"{soi_coh_path} doesn't exist ,please run create_soi.py first!")
    sys.exit(1)


# Process the data
soi_mm = soi * 1000
super_sboi = boi_mm.copy()
super_sboi[super_sboi==0]= np.nan
soi_mm[soi_mm==0]= np.nan
super_sboi[np.isnan(super_sboi)] = soi_mm[np.isnan(super_sboi)]

# Process and scale coherence data
super_sboi_coh = boi_coh.copy()
super_sboi_coh[super_sboi_coh==0] = np.nan
soi_coh[soi_coh==0]=np.nan
super_sboi_coh[np.isnan(super_sboi_coh)] = soi_coh[np.isnan(super_sboi_coh)]

# Scale coherence data to 0-255
super_sboi_coh[super_sboi_coh==0] = np.nan 
super_sboi_coh = np.clip(super_sboi_coh * 255, 0, 255)


# Export the result to a GeoTIFF
export_to_tiff(output_path, super_sboi, boi_path)
export_to_tiff(output_coh_path, super_sboi_coh, boi_path)

print(BLUE + "Super-sbovldiff successfully created!" + ENDC)