#!/bin/bash
#this will clip the SLCs ...

# note, for SAREZ:
# 005D_05199_131313 100A_05236_141313

# you may get hgt value as:
# source $LiCSARpath/lib/LiCSAR_bash_lib.sh
# framedir=`pwd`
# frame=`basename $framedir`
# hgt=$LiCSAR_public/`track_from_frame $frame`/$frame/metadata/$frame.geo.hgt.tif
# gdallocationinfo $hgt $lon $lat

# but better do it through python, i.e.:
# import xarray as xr; import rioxarray; import os
# frame='083D_12636_131313'
# hgt=os.path.join(os.environ['LiCSAR_public'], str(int(frame[:3])), frame, 'metadata', frame+'.geo.hgt.tif')
# a=rioxarray.open_rasterio(hgt)
# lon=-71.377; lat=-36.863; radius_km=25/2; radius_deg=radius_km/111; resol=0.00027;
# # a.sel(lon=(lon-radius_deg, lon+radius_deg), lat=(lat+radius_deg, lat-radius_deg))
# medhgt=float(a.sel(x=(lon-radius_deg, lon+radius_deg), y=(lat+radius_deg, lat-radius_deg), method='nearest').median())
# print('clip_slc.sh', 'outdir', str(lon-radius_deg), lon+radius_deg, lat-radius_deg, lat+radius_deg, medhgt, resol)


# and then generate hires geo (if it doesn't exist)
# cd volclip
# DEMDIR=$LiCSAR_proc/..../DEM
# geodir='geo'; mkdir $geodir
# masterslcdir='RSLC/'$master
# outres=0.00027027 # for 30 m, as 30/111000
# in python:
# python3 -c "from LiCSAR_lib.coreg_lib import geocode_dem; \
# geocode_dem('"$masterslcdir"', '"$geodir"', '"$DEMDIR"' , '.', '"$master"', "$outres")"

# however, we sometimes get the DEM with 'empty spaces' - so let's linear-interpolate and then nearest-neighbour interpolate the processed DEM
# so better check, using rasdt_pwr, and if see holes, just... check
# line 218 in coreg_lib and change geocode(lutfine,demseg,str(demwidth),hgtfile,str(width),str(length),'2','0',logfile) from '2' -> '1' (nearest neigh)
# or reinterpolate through python if needed - happened once, thus only as a comment:
#hgt=np.fromfile('GEOC/lookangles/20170311.geo.hgt', dtype=np.float32).byteswap()
#hgt=hgt.reshape(1257,1543)
#from lics_processing import *
#hgt=interpolate_nans(hgt, method='nearest')
#hgt.byteswap().tofile('GEOC/lookangles/20170311.geo.hgt')

export LiCSAR_subsets=$LiCSAR_procdir/subsets

# ok, now time to generate ifgs and unws
#echo "warning, the outfolder should be unique name (sorry for that, it is due to ifg generator) - so use e.g. VOLCID_008A etc."
if [ -z $7 ]; then echo "parameters are:";
echo "clip_slc.sh OUTFOLDER lon1 lon2 lat1 lat2 hei resolution [process_ifg] [init_only]"
echo "so e.g. clip_slc.sh CLIPPED -28.36 -27.3 38.49 38.8 600 0.00027 [1] [0]"
echo "where processing_ifg=1 means process data to ifgs after clipping"
echo "      and init_only=1 means only clip the core data for the frame subset"
exit;
fi

if [ ! -d RSLC ]; then echo "you need to be in the frame proc folder, i.e. the one with RSLC folder"; exit; fi
process_rslcs=1
process_ifgs=1
process_geomlis=1

if [ ! -z $8 ]; then
process_ifgs=$8;
fi

if [ ! -z $9 ]; then
 if [ $9 == 1 ]; then
	process_rslcs=0;
	process_ifgs=0;
	process_geomlis=0;
 fi
fi


source $LiCSARpath/lib/LiCSAR_bash_lib.sh

dizdir=`pwd`
frame=`basename $dizdir`
demdir=$LiCSAR_procdir/`track_from_frame $frame`/$frame/DEM
#framedir=$LiCSAR_procdir/`track_from_frame $frame`/$frame
if [ ! -d $demdir ]; then demdir=`pwd`/DEM; fi  # maybe it is there, locally?
if [ ! -d $demdir ]; then 
 # maybe running from ARC? try copy:
 scp -r -i /home/home02/earmla/.ssh/id_jasmin xfer1.jasmin.ac.uk:/gws/nopw/j04/nceo_geohazards_vol1/projects/LiCS/proc/current/`track_from_frame $frame`/$frame/DEM .
fi
if [ ! -d $demdir ]; then echo "something wrong with the frame "$frame". Fix this please"; exit; fi
dempar=$demdir/dem_crop.dem_par
if [ `ls geo/20??????.hgt 2>/dev/null | wc -l` == 0 ]; then
  echo "ERROR: the frame "$frame" has no hgt file in the geo folder - this needs fixing"; exit;
fi
master=`basename geo/20??????.hgt | cut -d '.' -f1`
tmpdir=$LiCSAR_temp/$frame/temp

slcpar=SLC/$master/$master.slc.par
if [ ! -f $slcpar ]; then echo "the folder is empty, or no mosaic of ref epoch exists - exiting"; exit; fi

#'clip_slc.sh 72.510 72.845 38.130 38.365 3934 0.00027'
#/gws/nopw/j04/nceo_geohazards_vol1/projects/LiCS/proc/current/subsets/SAREZ/005D
outdir=$1
if [ -d $outdir ]; then echo "STOP, this directory already exists. skipping for now, as this need some more thinking"; exit; fi
mkdir -p $outdir 2>/dev/null


lon1=$2
lon2=$3
lat1=$4
lat2=$5
hei=$6
resol=$7   # in degrees, so e.g. 0.00027 for 30 m
rgl=`echo $resol"*111000/2.3" | bc`
azl=`echo $resol"*111000/14" | bc`
resol_m=`python -c "print(round("$resol"*111111))"`

if [ $hei == 0 ]; then
  echo "warning, no height information provided. trying auto-extract but it might fail"
# echo "getting the avg height"
hei=`python3 -c "import xarray as xr; import rioxarray; import os; frame='"$frame"'; \
hgt=os.path.join(os.environ['LiCSAR_public'], str(int(frame[:3])), frame, 'metadata', frame+'.geo.hgt.tif'); \
a=rioxarray.open_rasterio(hgt); medhgt=float(a.sel(x=("$lon1","$lon2"), y=("$lat1", "$lat2"), method='nearest').median()); \
print(medhgt)"`
fi

# frame='083D_12636_131313'
# hgt=os.path.join(os.environ['LiCSAR_public'], str(int(frame[:3])), frame, 'metadata', frame+'.geo.hgt.tif')
# a=rioxarray.open_rasterio(hgt)
# lon=-71.377; lat=-36.863; radius_km=25/2; radius_deg=radius_km/111; resol=0.00027;
# # a.sel(lon=(lon-radius_deg, lon+radius_deg), lat=(lat+radius_deg, lat-radius_deg))
# medhgt=float(a.sel(x=(lon-radius_deg, lon+radius_deg), y=(lat+radius_deg, lat-radius_deg), method='nearest').median())





# e.g. for cz:
# clip_slc.sh czclip 18.51884 18.6357 49.7937 49.8515 293.784423828125 0.00027
rm corners_clip.tmp 2>/dev/null
cp $outdir/corners_clip.$frame corners_clip.tmp 2>/dev/null
mkdir -p $outdir/RSLC $outdir/log 2>/dev/null

if [ ! -f corners_clip.tmp ]; then
coord_to_sarpix $slcpar - $dempar $lat1 $lon1 $hei | grep "SLC/MLI range, azimuth pixel (int)" > corners_clip.tmp
coord_to_sarpix $slcpar - $dempar $lat2 $lon2 $hei | grep "SLC/MLI range, azimuth pixel (int)" >> corners_clip.tmp
coord_to_sarpix $slcpar - $dempar $lat1 $lon2 $hei | grep "SLC/MLI range, azimuth pixel (int)" >> corners_clip.tmp
coord_to_sarpix $slcpar - $dempar $lat2 $lon1 $hei | grep "SLC/MLI range, azimuth pixel (int)" >> corners_clip.tmp
fi

if [ ! -f $outdir/corners_clip.$frame ]; then
 cp corners_clip.tmp $outdir/corners_clip.$frame
fi

azi1=`cat corners_clip.tmp | rev | gawk {'print $1'} | rev | sort -n | head -n1`
azi2=`cat corners_clip.tmp | rev | gawk {'print $1'} | rev | sort -n | tail -n1`
let azidiff=azi2-azi1+1

rg1=`cat corners_clip.tmp | rev | gawk {'print $2'} | rev | sort -n | head -n1`
rg2=`cat corners_clip.tmp | rev | gawk {'print $2'} | rev | sort -n | tail -n1`
let rgdiff=rg2-rg1+1

if [ ! -f $outdir/RSLC/$master/$master.rslc ]; then
 mkdir -p $outdir/RSLC/$master
 # slc here is the master slc
 slc=`ls SLC/$master/$master.slc 2>/dev/null`
 if [ -z $slc ]; then 
	 mkdir -p $tmpdir
	 slc=$tmpdir/SLC/$master/$master.slc
	 mkdir -p $tmpdir/SLC/$master
	 tab=$tmpdir/SLC/$master/$master.tab
	 createSLCtab_frame SLC/$master/$master slc $frame > $tab
	 echo "mosaicing"
	 SLC_mosaic_ScanSAR $tab $slc $slc.par 20 4 0 >/dev/null
	 slcpar=$slc.par
 fi
 echo "clipping reference epoch"
 SLC_copy $slc $slcpar $outdir/RSLC/$master/$master.rslc $outdir/RSLC/$master/$master.rslc.par - - $rg1 $rgdiff $azi1 $azidiff - - >/dev/null 2>/dev/null
 x=$master
 multi_look $outdir/RSLC/$x/$x.rslc $outdir/RSLC/$x/$x.rslc.par $outdir/RSLC/$x/$x.rslc.mli $outdir/RSLC/$x/$x.rslc.mli.par $rgl $azl >/dev/null 2>/dev/null
fi

# now init
geodir=$outdir/geo.$resol_m'm'
if [ ! -d $geodir ]; then
	origdir=`pwd`
	echo "initializing the subset"
	cd $outdir
	if [ ! -d SLC/$master ]; then
	 mkdir -p SLC/$master
	 for x in `ls RSLC/$master/*`; do ln -s `pwd`/$x `pwd`/`echo $x | sed 's/RSLC/SLC/' | sed 's/rslc/slc/'`; done
	fi
	# prepare the geo folder

	#geodir='geo'; 
	mkdir -p $geodir
	masterslcdir='RSLC/'$master
	rm log/geo.err 2>/dev/null

	# in python:
	echo "preparing geocoding tables"
	python3 -c "from LiCSAR_lib.coreg_lib import geocode_dem; \
	 geocode_dem('"$masterslcdir"', '"$geodir"', '"$demdir"' , '.', '"$master"', "$resol")" > log/geo.log 2> log/geo.err

	if [ `grep -c 'Something' log/geo.err` -gt 0 ]; then 
		echo "some error in DEM fitting, skipping it now (might keep some slight DEM/geocoding shift)"
		rm -r $geodir; mkdir $geodir

		python3 -c "from LiCSAR_lib.coreg_lib import geocode_dem; \
		 geocode_dem('"$masterslcdir"', '"$geodir"', '"$demdir"' , '.', '"$master"', "$resol", skip_fit = True)"
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

	#echo $frame >> sourceframe.txt
	echo "clip_slc.sh "$outdir $lon1 $lon2 $lat1 $lat2 $hei $resol > sourcecmd.txt
	chmod -R 775 $outdir 2>/dev/null
	
	rm geo; rmdir GEOC #; mv GEOC GEOC.$resol_m'm'
	cd $origdir
	
	
fi

if [ $process_rslcs == 1 ]; then
	echo "performing full clipping"
	# ok, now clip the mosaics

	for x in `ls RSLC | grep 20`; do 
	 if [ -f RSLC/$x/$x.rslc ]; then
	 if [ ! -d $outdir/RSLC/$x ]; then
	   echo "clipping "$x
	   mkdir -p $outdir/RSLC/$x
	   SLC_copy RSLC/$x/$x.rslc RSLC/$x/$x.rslc.par $outdir/RSLC/$x/$x.rslc $outdir/RSLC/$x/$x.rslc.par - - $rg1 $rgdiff $azi1 $azidiff - - >/dev/null 2>/dev/null
	   # no need for multilooking here?... 
	   #multi_look $outdir/RSLC/$x/$x.rslc $outdir/RSLC/$x/$x.rslc.par $outdir/RSLC/$x/$x.rslc.mli $outdir/RSLC/$x/$x.rslc.mli.par $rgl $azl >/dev/null 2>/dev/null
	   # create_geoctiffs_to_pub.sh -M `pwd` $x >/dev/null   # to be improved
	 fi
	 fi
	done
fi

#procdir=$LiCSAR_procdir/`track_from_frame $frame`/$frame
#mkdir -p $procdir/subsets
#echo "now copy it to: "$LiCSAR_subsets
#finalout=$LiCSAR_subsets/$outdir
#echo "e.g. as: "$finalout
#echo "ln -s $finalout $procdir/subsets"





if [ $process_ifgs == 1 ]; then
	cd $outdir
	if [ ! -d geo ]; then ln -s geo.$resol_m'm' geo; fi
	if [ ! -d GEOC ]; then ln -s GEOC.$resol_m'm' GEOC; fi
	if [ ! -d GEOC.MLI ]; then ln -s GEOC.MLI.$resol_m'm' GEOC.MLI; fi
	# generate 'standard' connections ifgs
	echo "processing ifgs"
	framebatch_gapfill.sh -l -P 5 120 $rgl $azl
	echo "wait a bit and check.. tomorrow... for GEOC outputs"
fi

if [ $process_geomlis == 1 ]; then
	if [ ! -d geo ]; then ln -s geo.$resol_m'm' geo; fi
	if [ ! -d GEOC ]; then ln -s GEOC.$resol_m'm' GEOC; fi
	if [ ! -d GEOC.MLI ]; then ln -s GEOC.MLI.$resol_m'm' GEOC.MLI; fi
	echo 'generating MLI geotiffs'
	for x in `ls RSLC | grep 20`; do
	 if [ ! -d GEOC.MLI/$x/$x.geo.mli.tif ]; then
	   if [ ! -f $outdir/RSLC/$x/$x.rslc.mli ]; then
	     multi_look $outdir/RSLC/$x/$x.rslc $outdir/RSLC/$x/$x.rslc.par $outdir/RSLC/$x/$x.rslc.mli $outdir/RSLC/$x/$x.rslc.mli.par $rgl $azl >/dev/null 2>/dev/null
	   fi
	   #echo "geocoding "$x
	   create_geoctiffs_to_pub.sh -M `pwd` $x >/dev/null   # to be improved
	 fi
	done
fi

if [ ! -f $outdir/local_config.py ]; then
 echo "azlks="$azl > $outdir/local_config.py
 echo "rglks="$rgl >> $outdir/local_config.py
 echo "outres="$resol >> $outdir/local_config.py
 echo "resol_m="$resol_m >> $outdir/local_config.py
fi

chmod 777 $outdir

rm -r $tmpdir 2>/dev/null
cd $outdir; rm geo GEOC GEOC.MLI gmt.history 2>/dev/null; cd $dizdir
