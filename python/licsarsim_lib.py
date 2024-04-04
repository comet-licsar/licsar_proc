#!/usr/bin/env python3
# Simulation of SAR Intensity based on DEM
# M. Lazecky, 2024 (for Shailza Sharma, DEEPVOLC postdoc)

'''
Steps:
for one temporal epoch:
Inputs:
 - a DEM (geotiff) - will be simulated elsewhere
 - heading, centre inc angle, inc angle spread (e.g. +-5 deg) - target rg/azi sampling? (maybe not)
 - other params (corresponding to Sentinel-1 IW 'mid-swath')
Outputs:
 - geocoded simulated intensity

This will probably:
 - estimate satellite position (update SOVs)
 - convert DEM to radar coords
 - create geocoding tables for this DEM
 - simulate rdc intensity
 - geocode the simulated intensity
'''

# better to use py_gamma directly, but will just run direct commands
import os

def runcmd(cmd, message):
    print(message)
    os.system(cmd + ' >/dev/null 2>/dev/null')


