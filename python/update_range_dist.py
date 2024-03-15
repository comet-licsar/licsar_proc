#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Feb 12 13:10:34 2024

@author: cmagnard
"""

# import os
import glob
import py_gamma as pg

# os.chdir("/mydir")
par_list = glob.glob("*.par")
if not par_list:
    print('no par files found - are you in SLC/20xxxxxx directory?')
else:
    for par_file in par_list:
        par = pg.ParFile(par_file)
        r0 = par.get_value('near_range_slc', dtype = float, index = 0)
        rps = par.get_value('range_pixel_spacing', dtype = float, index = 0)
        nr = par.get_value('range_samples', dtype = int)
        r2 = r0 + (nr-1) * rps
        r1 = (r0 + r2)/2.0
        par.set_value('center_range_slc', r1, index = 0)
        par.set_value('far_range_slc', r2, index = 0)
        par.write_par(par_file)
        pg.update_par(par_file, par_file)
    print('done')