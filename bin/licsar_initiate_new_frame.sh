#!/bin/bash
#a script that will generate master, DEM  and all other necessary files to initiate a new frame
#module load licsar_proc #or _testing
curdir=$LiCSAR_procdir

dryrun=0; lastdays=365
setupmasterextra=''
outres=0.001
r=20
a=4
dolocal=0
tienshan=0
clip=0
sm=0 # stripmap frame

if [ -z $1 ];
then
 echo "Usage: licsar_initiate_new_frame.sh FRAME [MASTER_YYMMDD] "
 echo "where frame can be e.g. 010D_11058_131313 20180722"
 echo "(if master date is not given, it will choose automatically, from last 3 months data)"
 echo "(to include custom downloaded master files, don't forget to arch2DB.py them first)"
 echo "parameters:"
 echo " -H - do master in high resolution (r=5, a=1, res approx 15 m) - auto-applied if H is in framename"
 echo " -M - do master in medium resolution (towards 56 m outputs)"
 echo " -R 0.1 - custom resolution (here towards 0.1 deg) - note multilooking remains a=4,r=20, unless changed"
 echo " -a,-r - custom multilooking factors"
 echo " -D /path/to/dem.tif - use custom DEM - you may want to use gdal_merge.py -a_nodata -32768 .."
 echo " -V 365 would only output possible master epoch candidates..for last 365 days"
 echo " -L 365 would do same as -V but also auto-choose and process such master epoch"
 echo " -T - would include some extra Tien Shan related tuning"
 echo " -C lat1/lat2/lon1/lon2 - would establish a crop area for the frame"
 exit
fi

#improved getopts, finally
while getopts ":HMR:a:r:TD:V:L:C:" option; do
 case "${option}" in
  H) a=1; r=5; outres=0.00015; dolocal=1; echo "high resolution option enabled"
     ;;
  M) outres=0.0005; dolocal=1; echo "medium (56 m) resolution option enabled"
     ;;
  R) outres=$OPTARG; dolocal=1;
     ;;
  a) a=$OPTARG; dolocal=1;
     ;;
  r) r=$OPTARG; dolocal=1;
     ;;
  D) setupmasterextra=$setupmasterextra" -D "$OPTARG;
     #shift
     ;;
  T) tienshan=1; outres=0.0005; dolocal=1
     ;;
  V) dryrun=1;
     lastdays=$OPTARG;
     ;;
  L) lastdays=$OPTARG;
     setupmasterextra=$setupmasterextra" -L "$lastdays;
     ;;
  C) clip=1;
     cliparea=$OPTARG;
     dolocal=1;
     ;;
 esac
done
#shift
shift $((OPTIND -1))

if [ ! -z $2 ]; then
 getmaster="-m "$2
 else
 getmaster="-A 1"
fi


 frame=$1
 # update if the frame is for H:
 if [ ${frame:9:1} == 'H' ]; then
  a=1; r=5; outres=0.00015; dolocal=1;
 fi
 if [ `echo $frame | cut -d '_' -f 2` == 'SM' ]; then
   sm=1 #stripmap
   a=7; r=11; outres=0.0002777; dolocal=1;
   echo "Stripmap settings - overriding to"
   echo "a=7; r=11; outres=0.0002777"
 fi
 tr=`echo $frame | cut -d '_' -f1 | sed 's/^0//' | sed 's/^0//' | rev | cut -c 2- | rev`
 rmdir $curdir/$tr/$frame 2>/dev/null
 if [ -d $curdir/$tr/$frame ]; then
  echo "This frame already exists! Stopping here"
  echo "Check and remove(?) "$curdir/$tr/$frame
  exit
 fi
 if [ `echo $frame | grep -o '_' | wc -l` != 2 ]; then
  echo "Wrong frame name. Will try to continue (if not sure, cancel by CTRL-C)"
  #exit
  sleep 2
 fi

mkdir -p $curdir/$tr/$frame
cd $curdir/$tr/$frame


if [ $dryrun -gt 0 ]; then
 #lastdays=365
 echo "Finding master candidates in last "$lastdays" days."
 LiCSAR_setup_master.py -f $frame -d $curdir/$tr/$frame $getmaster -V $lastdays -r $r -a $a -o $outres $setupmasterextra
 exit
fi

if [ $dolocal == 1 ]; then
     echo "rglks = "$r > local_config.py
     echo "azlks = "$a >> local_config.py
     echo "outres = "$outres >> local_config.py
     if [ $clip == 1 ]; then
       echo "cliparea = "$cliparea >> local_config.py
     fi
fi
if [ $tienshan == 1 ]; then
    echo "tienshan = "1 >> local_config.py
fi

echo "Setting the master image and DEM for frame "$frame
LiCSAR_setup_master.py -f $frame -d $curdir/$tr/$frame $getmaster -r $r -a $a -o $outres $setupmasterextra
if [ ! -d $curdir/$tr/$frame/SLC ]; then
 echo "Something got wrong with the initiation"
 cd - 2>/dev/null
elif  [ -z `ls SLC` ]; then
 echo "Something got wrong with the initiation - no ref slc was created"
 cd - 2>/dev/null
else
 mkdir $curdir/$tr/$frame/LUT
 chmod 775 $curdir/$tr/$frame/LUT
 m=`ls SLC`
 if [ ! -d RSLC/$m ]; then 
  mkdir -p RSLC/$m
  for slcfile in `ls $curdir/$tr/$frame/SLC/$m/*`; do ln -s $slcfile `echo $slcfile | sed 's/SLC/RSLC/' | sed 's/slc/rslc/'`; done
 fi
 # 2024 - fix GAMMA error in centre_range
 cd SLC/$m; update_range_dist.py; cd -
 LiCSAR_05_mk_angles_master
 echo "Generating E-N-U files"
 if [ ! -f $curdir/$tr/$frame/geo/$m.lt_fine ]; then ln -s $curdir/$tr/$frame/geo/$m.lt $curdir/$tr/$frame/geo/$m.lt_fine; fi
 submit_lookangles.py -f $frame -t $tr
 echo "Generating land mask (would remove heights below 0 m)"
 landmask=$curdir/$tr/$frame/geo/landmask
 hgtgeo=$LiCSAR_public/$tr/$frame/metadata/$frame.geo.hgt.tif
 gmt grdlandmask -G$landmask.tif=gd:GTiff -R$hgtgeo -Df -N0/1/0/1/0
 #gmt grdconvert $landmask.grd $landmask.tif
 #rm $landmask.grd
 cp $landmask.tif $LiCSAR_public/$tr/$frame/metadata/$frame.geo.landmask.tif

#but this is wrong..so fixing.. ugly fast way:
# hgt=$LiCSAR_public/$tr/$frame/metadata/$frame.geo.hgt.tif
# lm=$LiCSAR_public/$tr/$frame/metadata/$frame.geo.landmask.tif
# g=`gmt grdinfo -C $hgt`
# ulx=`echo $g | gawk {'print $2'}`
# uly=`echo $g | gawk {'print $5'}`
# lrx=`echo $g | gawk {'print $3'}`
# lry=`echo $g | gawk {'print $4'}`
# #gdal_translate -of GTIFF -projwin $ulx $uly $lrx $lry $LiCSAR_procdir/GIS/GLDASp4_landmask_025d.nc4 $out
# ncpdq -O -h -a -lat $lm $lm.temp
# rm $lm
# gdal_translate -co COMPRESS=DEFLATE -co PREDICTOR=3 -a_ullr $ulx $uly $lrx $lry -of GTiff -a_srs EPSG:4326 $lm.temp $lm
# rm $lm.temp

 
 echo "Generating master MLI geotiff"
 create_geoctiffs_to_pub.sh -M `pwd` $m
 mkdir -p $LiCSAR_public/$tr/$frame/epochs 2>/dev/null
 mv GEOC.MLI/$m $LiCSAR_public/$tr/$frame/epochs/.
 rmdir GEOC.MLI
 
 echo "Generating public metadata file"
 submit_to_public_metadata.sh $frame
 #sometimes .xy is not generated..
 if [ ! -f $curdir/$tr/$frame/frame.xy ]; then cp $curdir/$tr/$frame/$frame'-poly.txt' $curdir/$tr/$frame/frame.xy; fi
 cp $curdir/$tr/$frame/$frame'-poly.txt' $LiCSAR_public/$tr/$frame/metadata/.
 echo "cleaning"
 rm -f $curdir/$tr/$frame/SLC/*/2*T*.I*sl* 2>/dev/null
 #removing also the mosaic 
 #rm -f $curdir/$tr/$frame/SLC/*/
 if [ $sm == 0 ]; then
   # do not remove for stripmaps
   rm -f $curdir/$tr/$frame/SLC/*/2???????.slc 2>/dev/null
 fi
 echo "done"

if [ $clip == 1 ]; then
 echo "clipping to requested area - WARNING, MUST BE lon1<lon2 etc"
 ulx=`echo $cliparea | cut -d '/' -f3`
 uly=`echo $cliparea | cut -d '/' -f2`
 lrx=`echo $cliparea | cut -d '/' -f4`
 lry=`echo $cliparea | cut -d '/' -f1`
 for tif in `ls $LiCSAR_public/$tr/$frame/metadata/*tif $LiCSAR_public/$tr/$frame/epochs/*/*tif`; do
   gdal_translate -projwin $ulx $uly $lrx $lry -co "COMPRESS=DEFLATE" -of GTiff -a_srs epsg:4326 $tif $tif.clip.tif
   mv $tif.clip.tif $tif
 done
 # just to clean the small png preview file (that were not clipped)
 rm $LiCSAR_public/$tr/$frame/epochs/*/*png
fi
fi

echo "checking for volcanoes in the frame (auto-init if the same relorb was not initialised yet)"
python3 -c "import volcdb; volcdb.initialise_subsets_in_frame('"$frame"')"

echo "changing permissions"
mkdir $LiCSAR_public/$tr/$frame/interferograms
chmod -R 775 $curdir/$tr/$frame $LiCSAR_public/$tr/$frame
chgrp -R gws_lics_admin $curdir/$tr/$frame
chgrp -R gws_lics_admin $LiCSAR_public/$tr/$frame
chmod 777 $LiCSAR_public/$tr/$frame/* # as some users outside the admin group will use store_to_curdir.sh
