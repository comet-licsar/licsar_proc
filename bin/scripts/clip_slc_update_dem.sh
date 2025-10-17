#!/bin/bash
# to change the DEM - included in clip_slc.sh.

# now we will replace the DEM to Copernicus DEM
# or... use path to your DEM (as geotiff) - needs to be over ellipsoid!
# source $LiCSARpath/lib/LiCSAR_bash_lib.sh

# please run this inside the subset directory, e.g.
# cd /gws/nopw/j04/nceo_geohazards_vol1/projects/LiCS/proc/current/subsets/firat/043A

# the most important part - get the lonlats
lonlats=`cat sourcecmd.txt | gawk {'print $3, $4, $5, $6'}`
if [ -e $lonlats 2>/dev/null ]; then echo "something is wrong - are you in the subset dir?"; exit; fi

if [ -d backup_oldDEM ]; then mv backup_oldDEM backup_oldDEM.old; fi

if [ -d DEM ]; then mkdir backup_oldDEM; mv DEM backup_oldDEM/.; fi

mkdir DEM
if [ ! -z $1 ]; then
  extdem=$1
  if [ ! -f $extdem ]; then
     echo "ERROR - the extra param must be existing geotiff of a DEM"
     exit
  fi
  LiCSAR_01_mk_crop_extDEM DEM/dem_crop 0 $extdem
  demname=`basename $extdem | cut -d '.' -f 1`
else
  # get the Copernicus DEM
  mk_copdem DEM/dem_crop $lonlats
  demname="CopernicusDEM_30m"
fi


# in python:
echo "preparing geocoding tables"
master=`ls SLC | head -n 1` # should be only one
masterslcdir=RSLC/$master

# update only resolution version that was set by local_config.py
source local_config.py # must have resol_m

# we may now want different resolution and multilooking than earlier, so....
mli=RSLC/$master/$master.rslc.mli.par
reprocmli=0
if [ ! -f $mli ]; then
  reprocmli=1
else
 if [ ! `get_value $mli azimuth_looks` == $azlks ]; then reprocmli=1; fi
 if [ ! `get_value $mli range_looks` == $rglks ]; then reprocmli=1; fi
 smli=SLC/$master/$master.slc.mli.par
 if [ ! `get_value $smli azimuth_looks` == $azlks ]; then reprocmli=1; fi
 if [ ! `get_value $smli range_looks` == $rglks ]; then reprocmli=1; fi
fi

if [ $reprocmli == 1 ]; then
  echo "regenerating ref epoch MLI"
  rm RSLC/$master/$master.rslc.mli RSLC/$master/$master.rslc.mli.par
  rm SLC/$master/$master.slc.mli SLC/$master/$master.slc.mli.par
  multi_look RSLC/$master/$master.rslc RSLC/$master/$master.rslc.par RSLC/$master/$master.rslc.mli RSLC/$master/$master.rslc.mli.par $rglks $azlks >/dev/null 2>/dev/null
  cp RSLC/$master/$master.rslc.mli SLC/$master/$master.slc.mli
  cp RSLC/$master/$master.rslc.mli.par SLC/$master/$master.slc.mli.par
fi

geodir=geo.$resol_m'm'  #`ls geo.*m -d | head -n 1`
mv $geodir GEOC.MLI.$resol_m'm' GEOC.meta.$resol_m'm' backup_oldDEM/. 2>/dev/null
mkdir $geodir
#resol_m=`echo $geodir | cut -d '.' -f 2 | cut -d 'm' -f 1`
resol=$outres #`cat sourcecmd.txt | rev | gawk {'print $1'} | rev`

mkdir -p $geodir
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


echo "DEM="$demname > metadata.txt
echo "clip_slc_update_dem.sh "$1 >> sourcecmd.txt

