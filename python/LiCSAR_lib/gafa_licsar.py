#!/usr/bin/env python

import h5py
import rasterio
import datetime as dt
from concurrent.futures import ProcessPoolExecutor, ThreadPoolExecutor
from datetime import datetime,timedelta,timezone
from pathlib import Path
import re
import numpy as np
import sys, os, subprocess, shutil
import glob
import shlex
import threading
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

AUX_DOWNLOAD_LOCK = threading.Lock()
AUTO_XFER_SERVER = object()
XFER_SERVERS = (
    "xfer-vm-01.jasmin.ac.uk",
    "xfer-vm-03.jasmin.ac.uk",
    "xfer-ph-01.jasmin.ac.uk",
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


def _raw_download_dir(frame_dir, acquisition_date=None):
    """Choose the shared RAW download directory, or the frame RAW directory."""
    licsar_slc = os.environ.get("LiCSAR_SLC")
    if not licsar_slc:
        raw_dir = frame_dir / "RAW"
        if acquisition_date is not None:
            raw_dir /= acquisition_date
        return raw_dir, False

    xfcpath = os.environ.get("XFCPATH")
    if xfcpath and Path(xfcpath).is_dir():
        return (Path(xfcpath) / "RAW").resolve(), True
    return (Path(licsar_slc).parent / "RAW").resolve(), True


def _find_xfer_server():
    """Return the first accessible transfer host, or None when none is reachable."""
    ssh_options = [
        "-A",
        "-o", "UserKnownHostsFile=/dev/null",
        "-o", "StrictHostKeyChecking=no",
    ]
    for server in XFER_SERVERS:
        try:
            result = subprocess.run(
                ["ssh", *ssh_options, server, "true"],
                stdout=subprocess.DEVNULL,
                stderr=subprocess.DEVNULL,
                timeout=10,
                check=False,
            )
        except (OSError, subprocess.TimeoutExpired):
            continue
        if result.returncode == 0:
            return server
    return None


def _select_monthly_epochs(epochs):
    """Select the available epoch closest to the 15th of each year-month."""
    epochs_by_month = {}
    for epoch in epochs:
        epoch_date = epoch if isinstance(epoch, dt.date) else dt.datetime.strptime(
            str(epoch), "%Y%m%d"
        ).date()
        epochs_by_month.setdefault((epoch_date.year, epoch_date.month), []).append(epoch_date)

    return [
        min(month_epochs, key=lambda date: (abs(date.day - 15), date))
        .strftime("%Y%m%d")
        for _, month_epochs in sorted(epochs_by_month.items())
    ]


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
    convert_to_gamma_output
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
    if convert_to_gamma_output:
        slc_fldname = frame_dir / 'SLC' / acquisition_date
        slc_fldname.mkdir(parents=True, exist_ok=True)
        gafa_fldname = next(date_dir.glob('GAF*'), None)
        if gafa_fldname is None:
            raise FileNotFoundError(f"No GAFA output directory found in {date_dir}")
        convert_to_gamma(gafa_fldname, slc_fldname)

def convert_to_gamma(gafa_out_dir, slc_out_dir):
    slc_out_dir = Path(slc_out_dir)
    subprocess.run(
        ["par_S1_GAFA_SLC", str(gafa_out_dir)],
        cwd=slc_out_dir,
        check=True,
    )

    output_pattern = re.compile(
        r"^(?P<date>\d{8})\.(?:vv|vh|hh|hv)\."
        r"(?P<swath>iw[123])\.slc(?P<suffix>\.par|\.tops_par)?$",
        re.IGNORECASE,
    )
    for output_file in slc_out_dir.iterdir():
        match = output_pattern.fullmatch(output_file.name)
        if not match or not output_file.is_file():
            continue
        suffix = match.group("suffix") or ""
        if suffix.lower() == ".tops_par":
            suffix = ".TOPS_par"
        renamed_file = slc_out_dir / (
            f"{match.group('date')}.{match.group('swath').upper()}.slc{suffix}"
        )
        if renamed_file.exists():
            raise FileExistsError(f"Cannot rename {output_file} to existing file {renamed_file}")
        output_file.rename(renamed_file)
    return True


def download_date(
    acquisition_date,
    frame,
    frame_dir,
    aux_dir = "/gws/ssde/j25a/nceo_geohazards/vol1/S1_AUX",
    download_aux = True,
    aux_types = ["INS", "CAL"],
    keep_aux_zip = False,
    download_orbits = True, # or link
    orbit_producttype = "POEORB",
    xfer_server = AUTO_XFER_SERVER,
    raw_download_dir = None,
):
    print("=" * 80)
    print(f"Processing frame {frame} for {acquisition_date}")
    print("=" * 80)
    date_dir = frame_dir / 'RAW' / acquisition_date
    date_dir.mkdir(parents=True, exist_ok=True)
    if raw_download_dir is None:
        raw_download_dir, shared_download_dir = _raw_download_dir(frame_dir, acquisition_date)
    else:
        shared_download_dir = Path(raw_download_dir).resolve() != date_dir.resolve()
    raw_download_dir = Path(raw_download_dir)
    raw_download_dir.mkdir(parents=True, exist_ok=True)
    if xfer_server is AUTO_XFER_SERVER:
        xfer_server = (
            _find_xfer_server()
            if shared_download_dir and os.environ.get("LiCSAR_SLC")
            else None
        )
    aux_dir = Path(aux_dir)
    aux_dir.mkdir(parents=True, exist_ok=True)
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
    if rawfiles is False or rawfiles is None or rawfiles.empty:
        raise RuntimeError(
            f"No RAW products found in CDSE for {frame}/{acquisition_date}."
        )

    print_section(f"RAW products found ({len(rawfiles)})")
    print_items(rawfiles["title"])

    # --------------------------------------------------
    # Get SLC products for the same frame/date
    # --------------------------------------------------
    start_dt = dt.datetime.strptime(acquisition_date, "%Y%m%d")
    slcfiles = lq.get_frame_files_date(frame, start_dt, True)
    if not slcfiles:
        print('WARNING - no SLC files in LiCSInfo DB for '+str(start_dt))
        # try from cdse:
        slcfiles = s1.get_images_for_frame(frame, acquisition_date, acquisition_date, prodType='SLC', asf = False)

    print_section(f"SLC coverage ({len(slcfiles)} SLC files)")
    covered_raw_products = identify_raw_covering_slcs(slcfiles, rawfiles['title'].values)

    if not covered_raw_products:
        raise RuntimeError(f"No RAW products cover the SLC files for {frame}/{acquisition_date}.")

    print_section(f"RAW products selected ({len(covered_raw_products)})")
    print_items(covered_raw_products)

    # --------------------------------------------------
    # Download Sentinel-1 AUX files needed by GAFA
    # --------------------------------------------------
    with AUX_DOWNLOAD_LOCK:
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
    print_section(f"RAW zip downloads -> {raw_download_dir}")
    for product in covered_raw_products:
        product = str(product).strip()

        if not product:
            continue

        if product.endswith(".zip"):
            zip_name = product
        else:
            zip_name = product + ".zip"

        out_file = raw_download_dir / zip_name

        if out_file.exists():
            print(f"  - exists: {zip_name}")
        else:
            print(f"  - downloading: {zip_name}")
            try:
                if xfer_server:
                    wget_cdse_path = shutil.which("wget_cdse")
                    if wget_cdse_path is None:
                        raise FileNotFoundError(
                            "Could not locate wget_cdse in the local PATH to invoke it "
                            "by absolute path on the transfer server."
                        )
                    wget_cdse_path = str(Path(wget_cdse_path).resolve())
                    remote_command = (
                        f"cd {shlex.quote(str(raw_download_dir))} && "
                        f"{shlex.quote(wget_cdse_path)} {shlex.quote(zip_name)}"
                    )
                    subprocess.run(
                        [
                            "ssh",
                            "-A",
                            "-o", "UserKnownHostsFile=/dev/null",
                            "-o", "StrictHostKeyChecking=no",
                            xfer_server,
                            remote_command,
                        ],
                        check=True,
                    )
                else:
                    subprocess.run(
                        ["wget_cdse", zip_name],
                        cwd=str(raw_download_dir),
                        check=True,
                    )
            except subprocess.CalledProcessError as e:
                print(f"FAILED: {zip_name}")
                print(f"Return code: {e.returncode}")
                continue

        if shared_download_dir:
            frame_raw_file = date_dir / zip_name
            if not os.path.lexists(frame_raw_file):
                if not out_file.is_file():
                    raise FileNotFoundError(
                        f"Downloaded RAW file is not accessible locally: {out_file}"
                    )
                os.symlink(out_file.resolve(), frame_raw_file)
                print(f"  - linked to frame RAW: {frame_raw_file}")


def process_frame_daterange(frame,
                            startdate = dt.date(2026,1,1), enddate = dt.date(2026,2,1),
                            normal_focusing = False, overwrite = False,
                            perform_download = False, perform_focusing = False,
                            convert_to_gamma_after_focusing = True, nproc = 1,
                            test_monthlies = False):
    ''' This will focus raw files for given frame, at given date range.
    We expect the processing to be performed from BATCH_CACHE_DIR where you should have RAW folder ready - or set it to download here?
    normal_focusing: if True, it will perform ordinary focusing (not extended bursts)
    perform_download: if True, it will autodownload RAW files and AUX files
    perform_focusing: if False, it will not perform focusing (useful if only download is needed)
    convert_to_gamma_after_focusing: doesn't do anything if perform_focusing is False
    nproc: number of acquisitions to focus in parallel
    test_monthlies: for Elgan's test, we will select only one epoch per each month
    '''
    # get epochs to process
    dates = s1.get_epochs_for_frame(frame, startdate, enddate)
    if test_monthlies:
        print('Selecting one epoch per each month')
        dates = _select_monthly_epochs(dates)

    print('\n Selected epochs are: ')
    print(dates)
    # get reference height
    metafile = os.path.join(os.environ['LiCSAR_public'],str(int(frame[:3])),frame,'metadata', 'metadata.txt')
    try:
        refhei = float(grep1line('avg_height',metafile).split('=')[1])
    except:
        print('WARNING: error loading avg height from metafile - setting ref hei to 100 m')
        refhei = 100.0
    # dates = ["20250324", "20250405"]
    batchdir = Path(os.environ["BATCH_CACHE_DIR"])
    # raw_root = Path(f"{batchdir}/ESBOI_02062026")
    frame_dir = batchdir / frame
    if not overwrite:
        # if exists in SLC already, do not download:
        dates2 = []
        for acquisition_date in dates:
            fldname = frame_dir / 'SLC' / acquisition_date
            if fldname.exists():
                print(acquisition_date+' exists in SLC folder. Skipping')
            else:
                dates2.append(acquisition_date)
        dates = dates2
    if perform_download:
        general_raw_dir, shared_download_dir = _raw_download_dir(frame_dir)
        xfer_server = _find_xfer_server() if shared_download_dir else None
        if xfer_server:
            print(f"Downloading RAW products via {xfer_server} to {general_raw_dir}")
        elif shared_download_dir:
            print(f"No transfer host available; downloading RAW products locally to {general_raw_dir}")
        download_args = [
            (
                acquisition_date,
                frame,
                frame_dir,
                _raw_download_dir(frame_dir, acquisition_date)[0],
            )
            for acquisition_date in dates
        ]
        if nproc == 1:
            for acquisition_date, download_frame, download_frame_dir, raw_download_dir in download_args:
                download_date(
                    acquisition_date,
                    download_frame,
                    download_frame_dir,
                    xfer_server=xfer_server,
                    raw_download_dir=raw_download_dir,
                )
        else:
            with ThreadPoolExecutor(max_workers=nproc) as executor:
                futures = [
                    executor.submit(
                        download_date,
                        acquisition_date,
                        download_frame,
                        download_frame_dir,
                        xfer_server=xfer_server,
                        raw_download_dir=raw_download_dir,
                    )
                    for acquisition_date, download_frame, download_frame_dir in download_args
                ]
                for future in futures:
                    future.result()
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
        if nproc == 1:
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
                    convert_to_gamma_after_focusing
                )
        else:
            with ProcessPoolExecutor(max_workers=nproc) as executor:
                futures = [
                    executor.submit(
                        focus_date,
                        acquisition_date,
                        frame,
                        frame_dir,
                        config,
                        custom_config,
                        extract_zips,
                        overwrite,
                        output_suffix,
                        convert_to_gamma_after_focusing,
                    )
                    for acquisition_date in dates
                ]
                for future in futures:
                    future.result()
        print("Focusing finished. Result is in " + str(frame_dir))
    print('done')
    print('now, you may want to do:')
    hgtif = os.path.join(os.environ['LiCSAR_public'], str(int(frame[:3])),frame,'metadata', frame+'.geo.hgt.tif')
    hgt = rasterio.open(hgtif)
    l,b,r,t=hgt.bounds
    print(f'mkdir DEM; mk_copdem DEM/dem_crop {l-0.2} {r+0.2} {b-0.2} {t+0.2}')
    print('for x in `ls SLC`; do convert2Xcomplex.sh SLC/$x SCOMPLEX; done')
    print('LiCSAR_02_coreg.py [-i?] -f '+frame+' -d . -m $epoch -R -E -k') # this would regenerate geo folder..
    print('framebatch_gapfill.sh -b 1')
    # this will not work due to cross normal SLC with extended SLC
    # print('licsar_make_frame.sh -n -f '+frame+' 0 0 2013-01-01 2013-01-02; mkdir $BATCH_CACHE_DIR/'+frame+'/log')
    # print('convert2Xcomplex.sh $BATCH_CACHE_DIR/'+frame+'/RSLC/*')
    # print('cd $BATCH_CACHE_DIR; framebatch_postproc_coreg.sh -f '+frame+' 1; rm -r '+frame+'/LUT')


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