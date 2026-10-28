#!/usr/bin/env python

import h5py
import datetime as dt
from datetime import datetime,timedelta,timezone
from pathlib import Path
import numpy as np
import sys, os, subprocess, shutil
import glob
import xarray as xr
import s1data as s1
from LiCSAR_misc import grep1line
import LiCSquery as lq
import s1data
import orbit_lib
import zipfile
import gafa
from gdar.base.raster import Collection, DataRaster
# GAFA functionality
from gafa.processor import focus_s1_tops_products, get_s1_raw_reader
from gafa.io.gafa.reader import reader as read_gafa
from download_s1_aux import (
    choose_aux_product,
    download_file,
    ensure_safe_dir,
    list_aux_products,
    parse_product_or_date,
)


def link_or_copy_file(src, dst_dir):
    src = Path(src)
    dst = dst_dir / src.name

    if dst.exists():
        print(f"  - exists: {dst.name}")
        return dst

    try:
        os.symlink(src, dst)
        print(f"  - linked: {dst.name}")
    except OSError:
        shutil.copy2(src, dst)
        print(f"  - copied: {dst.name}")

    return dst


def identify_raw_covering_slc(slcfile, rawfiles):
    ''' Helper function to identify which of the raw files cover given slc file

    Args:
        slcfile: slc file identifier, e.g.
        rawfiles: list of raw filenames, e.g. [

    Returns:
        raw file id that covers the given slcfile

    Note, useful functions:
    slcs = s1.get_images_for_frame('002A_16128_292828', s1.dt.date(2026,3,6), s1.dt.date(2026,3,7), prodType='SLC', asf = False)
    raws = s1.get_images_for_frame('002A_16128_292828', s1.dt.date(2026,3,6), s1.dt.date(2026,3,7), prodType='RAW', asf = False)
    '''
    t1slc, t2slc = get_dates_from_filename(slcfile)
    t1slc = t1slc + dt.timedelta(seconds=4)  # add tolerance
    t2slc = t2slc - dt.timedelta(seconds=4) # add tolerance
    for rawfile in rawfiles:
        t1, t2 = get_dates_from_filename(rawfile)
        if (t1 <= t1slc and t2 >= t2slc):
            return rawfile # full coverage
    print('None of the raw files cover this SLC')
    return None


def identify_raw_covering_slcs(slcfiles, rawfiles):
    """
    Identify the de-duplicated RAW products needed to cover the given SLC files.
    """
    raw_covers = []
    seen = set()
    for slcfile in slcfiles:
        raw_cover = identify_raw_covering_slc(slcfile, rawfiles)
        print(f"  - {Path(str(slcfile)).name}")
        print(f"    covered by: {raw_cover}")
        if raw_cover is None or raw_cover in seen:
            continue
        raw_covers.append(raw_cover)
        seen.add(raw_cover)
    return raw_covers


def get_dates_from_filename(filename):
    ''' Helper function to get t1 and t2 from given S1 filename

    Args:
        filename:  filename, e.g. 'S1A_IW_RAW__0SDV_20250324T233145_20250324T233217_058453_073AF4_530B'

    Returns:
        t1, t2:  dt.datetime
    '''
    t1, t2 = filename.split('_')[5:7]
    t1 = dt.datetime.strptime(t1, "%Y%m%dT%H%M%S")
    t2 = dt.datetime.strptime(t2, "%Y%m%dT%H%M%S")
    return t1, t2


''' Following functionality was developed by Muhammet Nergizci, an extremely talented PhD student. R.I.P. 2026
'''
### FOCUSING


def print_section(title):
    print()
    print(f"-- {title} --")


def print_items(items):
    for item in items:
        print(f"  - {item}")


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


def extract_safe_products(date_dir):
    """Extract RAW zip files in date_dir when the matching .SAFE folder is missing."""
    for zip_file in sorted(date_dir.glob("S1*_RAW_*.zip")):
        safe_dir = date_dir / zip_file.with_suffix(".SAFE").name
        if safe_dir.exists():
            print(f"  - exists: {safe_dir.name}")
            continue
        print(f"  - extracting: {zip_file.name}")
        with zipfile.ZipFile(zip_file) as zf:
            zf.extractall(date_dir)


def collect_raw_products(date_dir):
    """Return the SAFE products available for one acquisition date."""
    return sorted(date_dir.glob("S1*_RAW_*.SAFE"))


def focus_date(
    acquisition_date,
    frame,
    frame_dir,
    config,
    custom_config,
    extract_zips,
    overwrite,
    output_suffix,
):
    print("=" * 80)
    print(f"Focusing frame {frame} for {acquisition_date}")
    print("=" * 80)

    date_dir = frame_dir / 'RAW' / acquisition_date
    if not date_dir.exists():
        raise FileNotFoundError(f"RAW date directory does not exist: {date_dir}")

    if extract_zips:
        print_section(f"Extract RAW zip files -> {date_dir}")
        extract_safe_products(date_dir)

    raw_products = collect_raw_products(date_dir)
    if not raw_products:
        raise RuntimeError(f"No RAW .SAFE products found in {date_dir}")

    print_section(f"RAW products selected ({len(raw_products)})")
    for raw_product in raw_products:
        print(f"  - {raw_product.name}")

    print_section(f"GAFA output -> {date_dir}")
    focus_s1_tops_products(
        raw_products,
        config,
        date_dir,
        custom_config=custom_config,
        overwrite=overwrite,
        verbose=True,
        timing=True,
        output_suffix=output_suffix,
    )


def download_date(
    acquisition_date,
    frame,
    frame_dir,
    aux_dir = "/gws/ssde/j25a/nceo_geohazards/vol1/S1_AUX",
    download_aux = True,
    aux_types = ["INS", "CAL"],
    keep_aux_zip = False,
    download_orbits = True, # or link
    orbit_producttype = "POEORB"
):
    print("=" * 80)
    print(f"Processing frame {frame} for {acquisition_date}")
    print("=" * 80)
    date_dir = frame_dir / acquisition_date
    date_dir.mkdir(parents=True, exist_ok=True)
    if not os.path.exists(aux_dir):
        os.mkdir(aux_dir)
    # --------------------------------------------------
    # Get RAW product table for this acquisition date
    # --------------------------------------------------
    rawfiles = s1data.get_images_for_frame(
        frame,
        acquisition_date,
        acquisition_date,
        prodType="RAW",
        asf=False,
        outAspd=True
    )

    print_section(f"RAW products found ({len(rawfiles)})")
    print_items(rawfiles["title"])

    # --------------------------------------------------
    # Get SLC products for the same frame/date
    # --------------------------------------------------
    start_dt = dt.datetime.strptime(acquisition_date, "%Y%m%d")
    slcfiles = lq.get_frame_files_date(frame, start_dt, True)
    if not slcfiles:
        # try from cdse:
        slcfiles = s1.get_images_for_frame(frame, acquisition_date, acquisition_date, prodType='SLC', asf = False)

    print_section(f"SLC coverage ({len(slcfiles)} SLC files)")
    covered_raw_products = identify_raw_covering_slcs(slcfiles, rawfiles)

    if not covered_raw_products:
        raise RuntimeError(f"No RAW products cover the SLC files for {frame}/{acquisition_date}.")

    print_section(f"RAW products selected ({len(covered_raw_products)})")
    print_items(covered_raw_products)

    # --------------------------------------------------
    # Download Sentinel-1 AUX files needed by GAFA
    # --------------------------------------------------
    if download_aux:
        print_section("AUX files")
        aux_catalogs = {aux_type: list_aux_products(aux_type) for aux_type in aux_types}
        selected_aux = {}
        for product in covered_raw_products:
            product = str(product).strip()
            if not product:
                continue

            platform, acquisition_time = parse_product_or_date(product, None, None)
            for aux_type in aux_types:
                aux_product = choose_aux_product(
                    aux_catalogs[aux_type],
                    platform,
                    acquisition_time,
                    aux_type,
                )
                selected_aux[(aux_type, aux_product.name)] = aux_product

        for (aux_type, _), aux_product in sorted(selected_aux.items()):
            zip_file = aux_dir / f"{aux_product.name}.SAFE.zip"
            print(
                f"  - AUX_{aux_type}: {aux_product.name} "
                f"(valid from {aux_product.validity_start:%Y%m%dT%H%M%S})"
            )
            download_file(aux_product.download_url, zip_file)
            ensure_safe_dir(zip_file, aux_dir)
            if not keep_aux_zip:
                zip_file.unlink(missing_ok=True)

    # --------------------------------------------------
    # Link Sentinel-1 orbit files - needed by GAFA to reside in the aux_dir..
    # --------------------------------------------------
    if download_orbits:
        print_section("Orbit files")
        if orbit_lib is None:
            raise ImportError("download_orbits=True, but orbit_lib could not be imported.")

        selected_orbits = set()
        for product in covered_raw_products:
            product = str(product).strip()
            if not product:
                continue

            platform, acquisition_time = parse_product_or_date(product, None, None)
            acquisition_time = acquisition_time.replace(tzinfo=timezone.utc)
            orbit_list = orbit_lib.get_orbit_filenames_for_datetime(
                acquisition_time,
                producttype=orbit_producttype,
                s1ab=platform,
            )

            for orbit_file in orbit_list:
                selected_orbits.add(orbit_file)

        for orbit_file in sorted(selected_orbits):
            link_or_copy_file(orbit_file, aux_dir)

    # --------------------------------------------------
    # Download each RAW product into date_dir
    # --------------------------------------------------
    print_section(f"RAW zip downloads -> {date_dir}")
    for product in covered_raw_products:
        product = str(product).strip()

        if not product:
            continue

        if product.endswith(".zip"):
            zip_name = product
        else:
            zip_name = product + ".zip"

        out_file = date_dir / zip_name

        if out_file.exists():
            print(f"  - exists: {zip_name}")
            continue

        print(f"  - downloading: {zip_name}")
        try:
            subprocess.run(
                ["wget_cdse", zip_name],
                cwd=str(date_dir),
                check=True
            )
        except subprocess.CalledProcessError as e:
            print(f"FAILED: {zip_name}")
            print(f"Return code: {e.returncode}")
            continue


def process_frame_daterange(frame, startdate = dt.date(2026,1,1), enddate = dt.date(2026,2,1),
                            normal_focusing = False, overwrite = False,
                            perform_download = False, perform_focusing = False):
    ''' This will focus raw files for given frame, at given date range.
    We expect the processing to be performed from BATCH_CACHE_DIR where you should have RAW folder ready - or set it to download here?
    normal_focusing: if True, it will perform ordinary focusing (not extended bursts)
    perform_download: if True, it will autodownload RAW files and AUX files
    perform_focusing: if False, it will not perform focusing (useful if only download is needed)
    '''
    # get epochs to process
    dates = s1.get_epochs_for_frame(frame, startdate, enddate)
    print('\n Selected epochs are: ')
    print(dates)
    # get reference height
    metafile = os.path.join(os.environ['LiCSAR_public'],str(int(frame[:3])),frame,'metadata', 'metadata.txt')
    try:
        refhei = float(misc.grep1line('avg_height',metafile).split('=')[1])
    except:
        print('WARNING: error loading avg height from metafile - setting ref hei to 100 m')
        refhei = 100.0
    # dates = ["20250324", "20250405"]
    batchdir = Path(os.environ["BATCH_CACHE_DIR"])
    # raw_root = Path(f"{batchdir}/ESBOI_02062026")
    frame_dir = batchdir / frame
    if perform_download:
        for acquisition_date in dates:
            download_date(
                acquisition_date,
                frame,
                frame_dir
            )
        print('Downloading finished')
    if perform_focusing:
        # --------------------------------------------------
        # Config setup
        # --------------------------------------------------
        extract_zips = True
        if normal_focusing:
            config = os.path.join(os.environ['LiCSARpath'],'misc','gafa_s1_iw_ipf.toml')
            output_suffix = "_normal"
        else:
            config = os.path.join(os.environ['LiCSARpath'],'misc','gafa_s1_iw_fullspectrum.toml')
            output_suffix = "_extended"
        custom_config = {
            # "paths": {
            #    "sentinel1": {
            #        "aux": str(auxpath)
            #    }
            # },
            "geometry": {
                "ref_height": refhei
            }
        }
        for acquisition_date in dates:
            focus_date(
                acquisition_date,
                frame,
                frame_dir,
                config,
                custom_config,
                extract_zips,
                overwrite,
                output_suffix,
            )
        print("Focusing finished.")
    print('done')


'''

1.gafa_frame_inputs.py
2.gafa_focusing.py
are in:
/gws/ssde/j25a/nceo_geohazards/vol2/LiCS/temp/insar_proc/mnergizci/from_home/softwares/licsar_proc/python/LiCSAR_lib
but libs need:
/gws/ssde/j25a/nceo_geohazards/vol2/LiCS/temp/insar_proc/mnergizci/ESBOI_mynmar/examples/

'''

'''
to get slcs:
slcs = s1.get_images_for_frame('002A_16128_292828', s1.dt.date(2026,3,6), s1.dt.date(2026,3,7), prodType='SLC', asf = False)
slcs.sort()
raws = s1.get_images_for_frame('002A_16128_292828', s1.dt.date(2026,3,6), s1.dt.date(2026,3,7), prodType='RAW', asf = False)
raws.sort()
a=[]
for s in slcs:
    print(s)
    a.append(identify_raw_covering_slc(s, raws))

a.sort()
for x in a:
    print(x)


to get required PODs:
import orbit_lib as o
o.get_orbit_filenames_for_datetime(pd.Timestamp(2026,3,7,20,55,5, tz='UTC'))

to get required AUXs:
download_s1_aux.py --product S1D_IW_RAW__0SSH_20260307T205606_20260307T205638_001793_0028CF_7196 --platform S1D --dry-run

to get DEM:
mk_copdem_listonly 002A_16128_292828.geo.hgt.tif
'''