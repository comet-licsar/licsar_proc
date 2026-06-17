#!/usr/bin/env python3
'''
use merge_with_shifts( list_of_tifs, out_tif )
e.g.
list_of_tifs = glob.glob('*azi.geo.tif')
list_of_tifs.sort()
'''
import numpy as np
from osgeo import gdal

# the below means - set all nans to NODATA. due to the way gdal does things..
NODATA = -9999

def read_array(path):
    ds = gdal.Open(path)
    arr = ds.ReadAsArray().astype(float)
    gt = ds.GetGeoTransform()
    proj = ds.GetProjection()
    return arr, gt, proj

def write_array(path, arr, gt, proj, nodata=None):
    driver = gdal.GetDriverByName("GTiff")
    out = driver.Create(path, arr.shape[1], arr.shape[0], 1, gdal.GDT_Float32)
    out.SetGeoTransform(gt)
    out.SetProjection(proj)

    band = out.GetRasterBand(1)
    if nodata is not None:
        band.SetNoDataValue(nodata)

    band.WriteArray(arr)
    out.FlushCache()
    out = None


def compute_overlap(a1, gt1, a2, gt2):
    # Compute pixel overlap window
    x1_min, x1_max = gt1[0], gt1[0] + a1.shape[1] * gt1[1]
    y1_max, y1_min = gt1[3], gt1[3] + a1.shape[0] * gt1[5]

    x2_min, x2_max = gt2[0], gt2[0] + a2.shape[1] * gt2[1]
    y2_max, y2_min = gt2[3], gt2[3] + a2.shape[0] * gt2[5]

    # Overlap bounding box
    x_overlap_min = max(x1_min, x2_min)
    x_overlap_max = min(x1_max, x2_max)
    y_overlap_max = min(y1_max, y2_max)
    y_overlap_min = max(y1_min, y2_min)

    if x_overlap_min >= x_overlap_max or y_overlap_min >= y_overlap_max:
        return None  # no overlap

    # Convert to pixel coordinates
    def pix(gt, x, y):
        px = round((x - gt[0]) / gt[1])
        py = round((y - gt[3]) / gt[5])
        return px, py

    p1_min = pix(gt1, x_overlap_min, y_overlap_max)
    p1_max = pix(gt1, x_overlap_max, y_overlap_min)

    p2_min = pix(gt2, x_overlap_min, y_overlap_max)
    p2_max = pix(gt2, x_overlap_max, y_overlap_min)

    a1_ovl = a1[p1_min[1]:p1_max[1], p1_min[0]:p1_max[0]]
    a2_ovl = a2[p2_min[1]:p2_max[1], p2_min[0]:p2_max[0]]
    # Force same shape (fix off-by-one errors)
    min_rows = min(a1_ovl.shape[0], a2_ovl.shape[0])
    min_cols = min(a1_ovl.shape[1], a2_ovl.shape[1])
    a1_ovl = a1_ovl[:min_rows, :min_cols]
    a2_ovl = a2_ovl[:min_rows, :min_cols]

    return a1_ovl, a2_ovl



def merge_with_shifts(tif_list, out_tif):
    base_arr, base_gt, base_proj = read_array(tif_list[0])
    base_arr = np.where(np.isnan(base_arr), NODATA, base_arr)
    write_array("temp_mosaic.tif", base_arr, base_gt, base_proj, nodata=NODATA)

    for tif in tif_list[1:]:
        arr, gt, proj = read_array(tif)
        arr = np.where(np.isnan(arr), NODATA, arr)

        ovl = compute_overlap(base_arr, base_gt, arr, gt)
        if ovl is None:
            print(f"No overlap with {tif}, skipping shift")
            shift = 0
        else:
            a1_ovl, a2_ovl = ovl
            mask = (a1_ovl != NODATA) & (a2_ovl != NODATA)
            shift = np.nanmedian((a1_ovl[mask] - a2_ovl[mask]))
            print(f"Shifting {tif} by {shift}")

        arr_shifted = np.where(arr != NODATA, arr + shift, NODATA)

        write_array("temp_shifted.tif", arr_shifted, gt, proj, nodata=NODATA)

        warp_opts = gdal.WarpOptions(
            format="GTiff",
            creationOptions=["COMPRESS=LZW"],
            srcNodata=NODATA,
            dstNodata=NODATA
        )

        temp_path = "temp_mosaic.tif"
        #if not os.path.exists(temp_path):
        #    gdal.Warp(temp_path, ["temp_shifted.tif"], options=warp_opts)
        #else:
        gdal.Warp(temp_path, [temp_path, "temp_shifted.tif"], options=warp_opts)

        base_arr, base_gt, base_proj = read_array(temp_path)
        base_arr = np.where(np.isnan(base_arr), NODATA, base_arr)  # maybe not needed?

    write_array(out_tif, base_arr, base_gt, base_proj, nodata=NODATA)
    #os.system('gdal_merge.py -n nan -a_nodata nan -o '+out_path+' 'out_path+'.temp.tif '+tif_list[0]
    print(f"Merged output written to {out_path}")
