#!/bin/bash
# to change the DEM - can include this in clip_slc.sh if needed...
# now we will replace the DEM to Copernicus DEM (can be changed to any if you generate DEM/dem_crop.dem yourself


# please run this inside the subset directory, e.g.
# cd /gws/nopw/j04/nceo_geohazards_vol1/projects/LiCS/proc/current/subsets/firat/043A

# the most important part - get the lonlats
lonlats=`cat sourcecmd.txt | gawk {'print $3, $4, $5, $6'}`
if [ -e $lonlats 2>/dev/null ]; then echo "something is wrong - are you in the subset dir?"; exit; fi

if [ -d backup_oldDEM ]; then mv backup_oldDEM backup_oldDEM.old; fi
mkdir backup_oldDEM
if [ -d DEM ]; then mv DEM backup_oldDEM/.; fi

# get the Copernicus DEM
mkdir DEM
mk_copdem DEM/dem_crop $lonlats


# in python:
echo "preparing geocoding tables"
master=`ls SLC | head -n 1` # should be only one
masterslcdir=RSLC/$master

# update only resolution version that was set by local_config.py
source local_config.py # must have resol_m
geodir=geo.$resol_m'm'  #`ls geo.*m -d | head -n 1`
mv $geodir GEOC.MLI.$resol_m'm' GEOC.meta.$resol_m'm' backup_oldDEM/. 2>/dev/null
mkdir $geodir
#resol_m=`echo $geodir | cut -d '.' -f 2 | cut -d 'm' -f 1`
resol=$outres #`cat sourcecmd.txt | rev | gawk {'print $1'} | rev`

python3 -c "from LiCSAR_lib.coreg_lib import geocode_dem; \
	 geocode_dem('"$masterslcdir"', '"$geodir"', 'DEM' , '.', '"$master"', "$resol")" > log/geo_update.log 2> log/geo_update.err

if [ `grep -c 'Something' log/geo_update.err` -gt 0 ]; then 
	echo "some error in DEM fitting, skipping it now (might keep some slight DEM/geocoding shift)"
	rm -r $geodir; mkdir $geodir

	python3 -c "from LiCSAR_lib.coreg_lib import geocode_dem; \
	geocode_dem('"$masterslcdir"', '"$geodir"', 'DEM' , '.', '"$master"', "$resol", skip_fit = True)"
	cd $geodir
	ln -s $master.lt $master.lt_fine
	cd -
fi

# create mli file and geo geotiffs if needed
rmdir geo GEOC GEOC.MLI 2>/dev/null
rm geo GEOC GEOC.MLI 2>/dev/null
ln -s $geodir geo
echo "geocoding master mli and hgt"
create_geoctiffs_to_pub.sh -M `pwd` $master
#; ln -s 

echo "testing now - do ENUs"
LiCSAR_05_mk_angles_master
frame=`ls corners_clip* | head -n 1 | cut -d '.' -f 2`
submit_lookangles.py -f $frame -t `track_from_frame $frame` -l
create_geoctiff_lookangles.sh `pwd` $master #>/dev/null
mv GEOC.MLI GEOC.MLI.$resol_m'm'

for tif in `ls GEOC/lookangles/*tif`; do
 mv $tif `echo $tif | sed 's/'$master'/'$frame'/'`
done
echo "Generating land mask"
landmask=GEOC/lookangles/$frame.geo.landmask.tif
hgtgeo=`ls GEOC/lookangles/*.geo.hgt.tif | head -n 1`
gmt grdlandmask -G$landmask=gd:GTiff -R$hgtgeo -Df -N0/1/0/1/0

echo "testing now - store to GEOC.meta and deleting from GEOC - hope will not miss this"
mkdir GEOC.meta.$resol_m'm'
mv GEOC/lookangles/*.tif GEOC.meta.$resol_m'm'/.
rm -r GEOC/lookangles GEOC/geo 2>/dev/null


echo "DEM=CopernicusDEM_30m" > metadata.txt


