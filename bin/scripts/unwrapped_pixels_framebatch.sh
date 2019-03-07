#!/bin/bash

# GMT % unwrapped pixels for LiCSAR-generated GeoTIFF's:
\rm unwrapped_pixel_percent.list
touch unwrapped_pixel_percent.list
for dir in `ls GEOC/`
do
unw=`ls GEOC/$dir/*geo.unw.tif`
#echo making netcdf grid for $unw
#gdal_translate $unw -of GMT -a_nodata 99999 $unw".grd" 
grdclip -Sr0/NaN $unw -Gtemp.grd
perc=`grdinfo -M temp.grd | grep nodes | sed 's/(//g' | sed 's/%)//g' | awk '{print 100-$4}'`
echo ${perc}% of the pixels unwrapped
echo $dir $perc >> unwrapped_pixel_percent.list
\rm temp.grd
done
