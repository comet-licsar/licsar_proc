#!/usr/bin/env python3
'''
Created on 10/02/2024
@author: M. Nergizci, C. Magnard, M. Lazecky
University of Leeds
Gamma Remote Sensing
contact: mr.nergizci@gmail.com
'''
import argparse
import numpy as np
import os
import py_gamma as pg
import sys
import subprocess
import shutil
import time
import multiprocessing as mp
import LiCSAR_misc as misc
from modules_sw_mn import *  # functions saved here
import lics_unwrap as lu
import pickle

#Just for colored outputs
BLUE = '\033[94m'
ORANGE = '\033[38;5;208m'
ENDC = '\033[0m'  # ANSI code to end formatting


# Function to create empty `iw` files once
def create_empty_iw_files():
    master_tab_name = create_tab_file(master, framedir, frame, type='SLC')
    master_tab = np.array(framepath_tab(pg.read_tab(master_tab_name), framedir))

    nrows = master_tab.shape[0]
    nr, naz = [], []

    for i in range(nrows):
        SLC_par = pg.ParFile(master_tab[i][1])
        nr.append(SLC_par.get_value('range_samples', dtype=int, index=0))
        naz.append(SLC_par.get_value('azimuth_lines', dtype=int, index=0))

    for i in range(nrows):
        file_name = f'empty.iw{1 + i}.slc'
        full_path = os.path.join(temp_file, file_name)
        if not os.path.exists(full_path):
            pg.create_array(full_path, nr[i], naz[i], 5, 0.0, 0.0)
            print(f"{file_name} created.")


# Mosaicking and multilooking function to be executed in parallel
def mosaic_and_multilook(epoch):
    stdout_log_path = os.path.join(temp_file, f"{epoch}_soi_mosaic.out")
    stderr_log_path = os.path.join(temp_file, f"{epoch}_soi_mosaic.err")

    ## Redirect stdout and stderr to log files
    with open(stdout_log_path, "w") as out_file, open(stderr_log_path, "w") as err_file:
        sys.stdout = out_file
        sys.stderr = err_file

        try:
            # Tab files
            SLC1_tab_name = create_tab_file(epoch, framedir, frame, type='RSLC')
            master_tab_name = create_tab_file(master, framedir, frame, type='SLC')
            master_tab_name = create_tab_file(master, framedir, frame, type='RSLC')
            # Variable names
            SLC1_tab_mod1_name = os.path.join(tab_folder, epoch + '_mod1.SLC_tab')
            SLC1_tab_mod2_name = os.path.join(tab_folder, epoch + '_mod2.SLC_tab')
            SLC1_mod1_name = os.path.join(RSLC_folder, epoch, epoch + '_mod1.slc')
            SLC1_mod2_name = os.path.join(RSLC_folder, epoch, epoch + '_mod2.slc')
            mli1_mod1_name = os.path.join(RSLC_folder, epoch, epoch + '_mod1.mli')
            mli1_mod2_name = os.path.join(RSLC_folder, epoch, epoch + '_mod2.mli')

            # Number of looks
            rlks, azlks = 20, 4

            # Read SLC_tab files
            SLC1_tab = np.array(framepath_tab(pg.read_tab(SLC1_tab_name), framedir))
            master_tab = np.array(framepath_tab(pg.read_tab(master_tab_name), framedir))

            # Create the full path for master_tab
            master_tab_long_name = os.path.join(framedir, 'tab', master + '_tab_long')
            if not os.path.exists(master_tab_long_name):
                with open(master_tab_long_name, 'w') as file:
                    for line in master_tab:
                        file.write(" ".join(line) + "\n")

            # mode1, mode2 creation
            SLC1_tab_mod1, SLC1_tab_mod2 = SLC1_tab.copy(), SLC1_tab.copy()
            SLC1_tab_mod1[1][0] = os.path.join(temp_file, 'empty.iw2.slc')
            SLC1_tab_mod2[0][0] = os.path.join(temp_file, 'empty.iw1.slc')
            SLC1_tab_mod2[2][0] = os.path.join(temp_file, 'empty.iw3.slc')

            # Save the new tabs to temp_file
            pg.write_tab(SLC1_tab_mod1, SLC1_tab_mod1_name)
            pg.write_tab(SLC1_tab_mod2, SLC1_tab_mod2_name)

            # Update burst_win parameters using ext_burst_win values
            nrows = master_tab.shape[0]
            for i in range(nrows):
                TOPS_par_master = pg.ParFile(master_tab[i][2])
                TOPS_par = pg.ParFile(SLC1_tab[i][2])
                nburst = TOPS_par.get_value('number_of_bursts', dtype=int, index=0)
                for b in range(nburst):
                    ext_burst_win = TOPS_par_master.get_value(f'ext_burst_win_{b + 1}')
                    TOPS_par.set_value(f'burst_win_{b + 1}', ext_burst_win[0], index=0)
                    TOPS_par.set_value(f'burst_win_{b + 1}', ext_burst_win[1], index=1)
                    TOPS_par.set_value(f'burst_win_{b + 1}', ext_burst_win[4], index=4)
                    TOPS_par.set_value(f'burst_win_{b + 1}', ext_burst_win[5], index=5)
                TOPS_par.write_par(SLC1_tab[i][2])
                pg.update_par(SLC1_tab[i][2], SLC1_tab[i][2])

            # Run mosaicking and multilooking
            if not os.path.exists(SLC1_mod1_name) or os.path.getsize(SLC1_mod1_name) == 0:
                pg.SLC_mosaic_ScanSAR(SLC1_tab_mod1_name, SLC1_mod1_name, SLC1_mod1_name + '.par', rlks, azlks, 0, master_tab_long_name)
            if not os.path.exists(SLC1_mod2_name) or os.path.getsize(SLC1_mod2_name) == 0:
                pg.SLC_mosaic_ScanSAR(SLC1_tab_mod2_name, SLC1_mod2_name, SLC1_mod2_name + '.par', rlks, azlks, 0, master_tab_long_name)

            if not os.path.exists(mli1_mod1_name) or os.path.getsize(mli1_mod1_name) == 0:
                pg.multi_look(SLC1_mod1_name, SLC1_mod1_name + '.par', mli1_mod1_name, mli1_mod1_name + '.par', rlks, azlks)
            if not os.path.exists(mli1_mod2_name) or os.path.getsize(mli1_mod2_name) == 0:
                pg.multi_look(SLC1_mod2_name, SLC1_mod2_name + '.par', mli1_mod2_name, mli1_mod2_name + '.par', rlks, azlks)

        except Exception as e:
            print(f"Error processing epoch {epoch}: {e}")

        finally:
            sys.stdout = sys.__stdout__
            sys.stderr = sys.__stderr__


### Main processing
if __name__ == "__main__":
    
    # Argument parser
    parser = argparse.ArgumentParser(description="Process mosaic and multilooking in parallel.")
    parser.add_argument("--n_para", type=int, default=10, help="Number of parallel processes (default: 10)")
    parser.add_argument("--epoch_list", type=str, help="Path to a file with a list of epochs to process (optional)")
    args = parser.parse_args()

    # Start time
    start_time = time.time()

    # Variables
    tempdir = os.getcwd()
    frame = os.path.basename(tempdir)
    # batchdir = os.environ['BATCH_CACHE_DIR']
    framedir = os.path.join(tempdir)
    tr = int(frame[:3])
    metafile = os.path.join(os.environ['LiCSAR_public'], str(tr), frame, 'metadata', 'metadata.txt')
    master = misc.grep1line('master=', metafile).split('=')[1]
    RSLC_folder = os.path.join(framedir, 'RSLC')
    tab_folder = os.path.join(framedir, 'tab')
    n_para = args.n_para
    ## Output directory
    temp_file = os.path.join(framedir, 'temp_data')
    if not os.path.exists(temp_file):
        os.makedirs(temp_file)

    
    # Create empty `iw` files once
    create_empty_iw_files()
    print('empty subswaths created...')
    # Load epochs from the list file or from RSLC_folder
    if args.epoch_list:
        with open(args.epoch_list, "r") as f:
            epochs = [line.strip() for line in f if line.strip()]
    else:
        epochs = [epoch for epoch in os.listdir(RSLC_folder) if os.path.isdir(os.path.join(RSLC_folder, epoch))]
        # if master in epochs: ## Keep it, I am not sure master epoch should be excluded or not.
        #     epochs.remove(master)
            
    # for epoch in epochs: ## MN this is for sequential processing for testing via breakpoint()
    #     print(f"{BLUE}Processing epoch: {epoch}{ENDC}")
    #     mosaic_and_multilook(epoch)

    # Mosaicking and multilooking in parallel 
    with mp.Pool(processes=n_para) as pool:
        pool.map(mosaic_and_multilook, epochs)

    end_time = time.time()
    execution_time = end_time - start_time
    print(f"Mosaicking and multilooking took {round(execution_time)} seconds!")