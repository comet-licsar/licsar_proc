#!/usr/bin/env python
# https://github.com/drf5n/drf5n-public/blob/master/gdalwarp2match.py


from osgeo import gdal, gdalconst
import argparse

parser = argparse.ArgumentParser(description='Use GDAL to reproject a raster to match the extents and res of a template')
parser.add_argument("source", help="Source file")
parser.add_argument("template", help = "template with extents and resolution to match")
parser.add_argument("destination", help = "destination file (geoTIFF)")
args = parser.parse_args()
print(args)



# Source
src_filename = args.source
src = gdal.Open(src_filename, gdalconst.GA_ReadOnly)
src_proj = src.GetProjection()
src_geotrans = src.GetGeoTransform()

# We want a section of source that matches this:
match_filename = args.template
match_ds = gdal.Open(match_filename, gdalconst.GA_ReadOnly)
match_proj = match_ds.GetProjection()
match_geotrans = match_ds.GetGeoTransform()
wide = match_ds.RasterXSize
high = match_ds.RasterYSize

# Output / destination
dst_filename = args.destination
dst = gdal.GetDriverByName('GTiff').Create(dst_filename, wide, high, 1, gdalconst.GDT_Float32)
dst.SetGeoTransform( match_geotrans )
dst.SetProjection( match_proj)

# Do the work
# GRA_Bilinear Bilinear give diamond-shaped artifacts
# GRA_CubicSpline Cubic-Spline smooths them
# GRA_NearestNeighbour nearest neighbor???? 
interpolation=gdalconst.GRA_NearestNeighbour
#gdal.ReprojectImage(src, dst, src_proj, match_proj, gdalconst.GRA_Bilinear)
#gdal.ReprojectImage(src, dst, src_proj, match_proj, gdalconst.GRA_CubicSpline)
gdal.ReprojectImage(src, dst, src_proj, match_proj, interpolation)

del dst # Flush

