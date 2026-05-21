#!/usr/bin/env python3
import sys
import os
import subprocess
import shutil
from datetime import timezone
from pathlib import Path

repo_root = Path(__file__).resolve().parent.parent
sys.path.append(str(repo_root))

# --------------------------------------------------
# LiCSAR python path and environment
# --------------------------------------------------
sys.path.append("/home/users/mnergiz/softwares/licsar_proc/python/LiCSAR_lib")
os.environ["LiCSAR_SLC"] = "/gws/ssde/j25a/nceo_geohazards/vol2/LiCS/temp/SLC/S1"
import s1data
import orbit_lib

from tools.download_s1_aux import (
    choose_aux_product,
    download_file,
    ensure_safe_dir,
    list_aux_products,
    parse_product_or_date,
)

# --------------------------------------------------
# Inputs
# --------------------------------------------------
frame = "106D_06844_131313"
startdate = "20250324"
enddate = "20250324"
batchdir= os.environ['BATCH_CACHE_DIR']
# raws_mynmar = Path(f"{batchdir}/RAW")
raws_mynmar = Path(f"{batchdir}/RAWS_mynmar")
out_dir = raws_mynmar / frame
# aux_dir = Path("/gws/ssde/j25a/nceo_geohazards/vol2/LiCS/temp/insar_proc/mnergizci/ESBOI_mynmar/S1_AUX")
aux_dir = Path("/gws/ssde/j25a/nceo_geohazards/vol1/S1_AUX") #MN: Milan please check this location, is there ok to save? gafa wants the orbit files in same folder so I just soft linked. 

# wget_cdse = Path("/home/users/mnergiz/softwares/licsar_proc/bin/scripts/wget_cdse")
download_aux = True
aux_types = ["INS", "CAL"] #PP1
keep_aux_zip = False
download_orbits = True # or link
orbit_producttype = "POEORB"

out_dir.mkdir(parents=True, exist_ok=True)
aux_dir.mkdir(parents=True, exist_ok=True)


def link_or_copy_file(src, dst_dir):
    src = Path(src)
    dst = dst_dir / src.name

    if dst.exists():
        print(f"Already exists, skipping: {dst.name}")
        return dst

    try:
        os.symlink(src, dst)
        print(f"Linked orbit: {dst}")
    except OSError:
        shutil.copy2(src, dst)
        print(f"Copied orbit: {dst}")

    return dst

# --------------------------------------------------
# Get RAW product table
# --------------------------------------------------
files = s1data.get_images_for_frame(
    frame,
    str(startdate),
    str(enddate),
    prodType="RAW",
    asf=False,
    outAspd=True
)

print(files[["title"]])

# --------------------------------------------------
# Download Sentinel-1 AUX files needed by GAFA
# --------------------------------------------------
if download_aux:
    aux_catalogs = {aux_type: list_aux_products(aux_type) for aux_type in aux_types}
    selected_aux = {}
    for product in files["title"]:
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
            f"Downloading AUX_{aux_type}: {aux_product.name} "
            f"(valid from {aux_product.validity_start:%Y%m%dT%H%M%S})"
        )
        download_file(aux_product.download_url, zip_file)
        ensure_safe_dir(zip_file, aux_dir)
        if not keep_aux_zip:
            zip_file.unlink(missing_ok=True)

# --------------------------------------------------
# Link Sentinel-1 orbit files needed by GAFA
# --------------------------------------------------
if download_orbits:
    if orbit_lib is None:
        raise ImportError("download_orbits=True, but orbit_lib could not be imported.")

    selected_orbits = set()
    for product in files["title"]:
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
# Download each RAW product into raws_mynmar/frame
# --------------------------------------------------
for product in files["title"]:
    product = str(product).strip()

    if not product:
        continue

    if product.endswith(".zip"):
        zip_name = product
    else:
        zip_name = product + ".zip"

    out_file = out_dir / zip_name

    if out_file.exists():
        print(f"Already exists, skipping: {zip_name}")
        continue

    print(f"Downloading: {zip_name}")
    try:
        subprocess.run(
            ["wget_cdse", zip_name],
            cwd=str(out_dir),
            check=True
        )
    except subprocess.CalledProcessError as e:
        print(f"FAILED: {zip_name}")
        print(f"Return code: {e.returncode}")
        continue

print("Finished.")
