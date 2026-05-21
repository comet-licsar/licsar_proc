#!/usr/bin/env python3

from pathlib import Path
try:
    import gafa  # noqa: F401
except:  # noqa: E722
    # This happens if the PYTHONPATH does not contain the GAFA module.
    import sys
    _gafa_path = Path().cwd().parent / "src"
    sys.path.append(_gafa_path.as_posix())
    import gafa  # noqa: F401

import numpy as np
from matplotlib import pyplot as plt

from gdar.base.raster import Collection, DataRaster

# GAFA functionality
from gafa.processor import focus_s1_tops_products, get_s1_raw_reader
from gafa.io.gafa.reader import reader as read_gafa
from gafa.rastertools import srd_mosaic
import os 

def read_gafa_slc(path) -> DataRaster:
    """
    Function that provides a reader for a folder of focused bursts in NetCDF format.
    """
    match path.stem[4:6]:
        case "IW":
            swaths = ["IW1", "IW2", "IW3"] 
        case "EW":
            swaths = ["EW1","EW2","EW3","EW4","EW5"]
    pols = [path.stem[15]]
    r = {}
    for sw in swaths:
        s = {}
        for p in pols:
            fslc_gafa = sorted([f for f in path.glob(f"**/*{sw}-{p}*.nc")])
            b = []
            for f in fslc_gafa:
                b.append(read_gafa(f))
            if len(b) != 0:
                s[p] = b
        if len(s) != 0:
            r[sw] = s
    return Collection(r) if len(r) > 0 else None

def multilook(d, my: int, mx: int, y0: int = 0, x0: int = 0):
    """
    Boxcar multilook filter with factors (my, mx).
    """
    ny, nx = d.shape[0] - y0, d.shape[1] - x0
    ly, lx = ny // my, nx // mx
    ky, kx = ly * my, lx * mx
    return np.mean(
        np.mean(d[y0 : y0 + ky, x0 : x0 + kx].reshape(ly, my, lx, mx), axis=-1), axis=-2
    )


# --------------------------------------------------
# input Path setup
# --------------------------------------------------
batchdir= os.environ['BATCH_CACHE_DIR']
rawpath = Path(f"{batchdir}/RAWS_mynmar")
auxpath = Path("/gws/ssde/j25a/nceo_geohazards/vol1/S1_AUX")
normal= True #True: focusing normal, False: focusing extended
# --------------------------------------------------
# Config setup
# --------------------------------------------------
if normal:
    config = Path.cwd().parents[0] / "gafa_config" / "s1" / "s1_iw_ipf.toml"
    opath= Path(f"{batchdir}/SLCS_normal")
else:
    config = Path.cwd().parents[0] / "gafa_config" / "s1" / "s1_iw_fullspectrum.toml"
    opath= Path(f"{batchdir}/SLCS_extended")

custom_config = {
    "paths": {
        "sentinel1": {
            "aux": str(auxpath)
        }
    },
    "geometry": {
        "ref_height": 100.0
    }
}

# --------------------------------------------------
# output Path setup
# --------------------------------------------------
opath.mkdir(parents=True, exist_ok=True)
frame="106D_06844_131313"

for raw in os.listdir(rawpath / frame):
    if raw.endswith(".SAFE"):
        print(raw)
        fraw = rawpath / frame / raw
        print(f"Processing {fraw}...")
        # # Run the focusing
        focus_s1_tops_products(fraw, config, opath, custom_config=custom_config, overwrite=True, verbose=True, timing=True)