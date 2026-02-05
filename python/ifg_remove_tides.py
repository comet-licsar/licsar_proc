#!/usr/bin/env python

import sys
import rioxarray as r
import numpy as np

# Read command-line arguments
hgtfile = sys.argv[1]  # Height file
ifile = sys.argv[2]    # Input interferogram
t1f = sys.argv[3]      # Time 1 file
t2f = sys.argv[4]      # Time 2 file
ofile = sys.argv[5]    # Output file
ocfile = sys.argv[6] if len(sys.argv) > 6 else None # Output corection file, if exists

# Load raster files
h = r.open_rasterio(hgtfile).squeeze()  # Remove unnecessary dimensions
i = r.open_rasterio(ifile).squeeze()
t1 = r.open_rasterio(t1f).squeeze()
t2 = r.open_rasterio(t2f).squeeze()

# Replace height values with interferogram values
# h = h.copy()
h.values = i.values
h=h.where(h!=0)

dt = h.copy()
# Compute dt based on file naming
if "bovl" in ifile:  # Check if 'bovl' is in filename
    dt.values = (t1.values - t2.values) * 1000 # Convert m2mm
    if ocfile:  # Only save if ocfile is provided
        dt.rio.to_raster(ocfile)
else:
    dt.values = (t1.values - t2.values) * 226.56 #Convert m2rad

# Apply phase correction
h.values = h.values - dt.values
if "bovl" not in ifile:
    h.values = np.angle(np.exp(1j * h.values))  # Convert phase to [-π, π]

# Save corrected raster
h.rio.to_raster(ofile)
