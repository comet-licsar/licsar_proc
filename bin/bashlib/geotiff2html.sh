#!/bin/bash

in=input.tif
outdir=AHB_VE
gmtcpt=
gmtcpt2gdalcpt.py $gmtcpt $gmtcpt.gdal
gdaldem color-relief $in $gmtcpt.gdal $in.coloured.tif
gdal2tiles.py $in.coloured.tif $outdir

# anebo:
cpt=SCM.roma_r
cpt=RdBu_r
in=VEL_asc_eur.tif
out=AHB_asc
LiCSBAS_color_geotiff.py -i $in -c $cpt --cmin -15 --cmax 15 -o $in.coloured.tif
LiCSBAS_color_geotiff2tiles.py -i $in.coloured.tif -o $out --n_para 4
