#!/bin/bash

source $LiCSARpath/lib/LiCSAR_bash_lib.sh
LOTUS=1

if [ -z $2 ]; then
  cat << End_of_Usage

  Purpose: Merge and unwrap interferogram pairs shared by two adjacent LiCSAR frames.

  Usage: frame_merge_unwrap.sh [-i IFG] [-l IFG.list] FRAME_ID1 FRAME_ID2

         e.g.
         frame_merge_unwrap.sh 106D_05447_131313 106D_05646_131313
         or
         frame_merge_unwrap.sh -i 20210911_20210929 106D_05447_131313 106D_05646_131313
         or
         frame_merge_unwrap.sh -l ifg.list 106D_05447_131313 106D_05646_131313

         A new frame ID will be given (e.g. 106D_05447_05646).

  Author: Jin Fang @ Uni Leeds.
	  Thanks to Milan for optimisation!
  Date: Nov, 2021

End_of_Usage
  exit
fi

FMU=1

while getopts 'i:l:' OPT; do
    case $OPT in
        i) FMU=0; echo $OPTARG > ifg_merge_unwrap.list
           ;;
        l) FMU=0; cat $OPTARG > ifg_merge_unwrap.list
           ;;
    esac
done
shift $((OPTIND -1))

frame1=$1
frame2=$2
track=`track_from_frame $frame1`
frame=${frame1:0:11}${frame2:5:5}

if [ -d $frame/interferograms ]; then
 if [ -d $frame/GEOC ]; then
  echo "check the merged frame folder - seems you have both old and new merge done?"
  exit
 fi
 mv $frame/interferograms $frame/GEOC
 mv $frame/metadata/* $frame/GEOC/. 2>/dev/null
 rmdir $frame/metadata 2>/dev/null
else
 mkdir -p $frame/GEOC
fi

if [ $FMU -eq 0 ]; then
   mv ifg_merge_unwrap.list $frame/
elif [ $FMU -eq 1 ]; then
   if [ -f $frame/ifg_merge_unwrap.list ]; then
       rm $frame/ifg_merge_unwrap.list
   fi
   echo "Checking the total number of interferogram pairs to process..."
   if [ ! -d $frame1/GEOC ]; then ln -s $frame1/interferograms $frame1/GEOC; fi
   if [ ! -d $frame1/GEOC ]; then ln -s $frame2/interferograms $frame2/GEOC; fi
   for i in `ls -d $frame1/GEOC/20*`;
   do
       if [ -d "$frame2/GEOC/${i:0-17:17}" ]; then
           ifg=${i:0-17:17}
           if [ ! -d $frame/GEOC/$ifg ]; then
             echo $ifg >> $frame/ifg_merge_unwrap.list
           fi
       fi
   done
fi

start_time=`date +%s`

sum=`cat $frame/ifg_merge_unwrap.list | wc -l`
echo "The number of interferogram pairs to be merged and unwrapped: $sum"

if [ -f $frame/ifg_failed.list ]; then
   rm $frame/ifg_failed.list
fi

echo "Preparing frame files"
mkdir -p $frame/GEOC
mask1=$LiCSAR_public/$track/$frame1/metadata/$frame1.geo.landmask.tif
mask2=$LiCSAR_public/$track/$frame2/metadata/$frame2.geo.landmask.tif


for more in E N U hgt landmask; do
if [ ! -f $frame/GEOC/$frame.geo.$more.tif ]; then
 echo "merging "$more
 morefrom1=$LiCSAR_public/$track/$frame1/metadata/$frame1.geo.$more.tif
 morefrom2=$LiCSAR_public/$track/$frame2/metadata/$frame2.geo.$more.tif
 gdal_merge.py -n 0 -co COMPRESS=DEFLATE -o $frame.geo.$more.tif $morefrom1 $morefrom2 # output file name containing path got error
 mv $frame.geo.$more.tif $frame/GEOC/.
fi
done

if [ ! -f $frame/GEOC/$frame.geo.mli.tif ]; then
m1=`get_master $frame1`
morefrom1=$LiCSAR_public/$track/$frame1/$m1/$m1.geo.mli.tif
m2=`get_master $frame2`
morefrom2=$LiCSAR_public/$track/$frame2/$m2/$m2.geo.mli.tif
if [ -f $morefrom1 ] && [ -f $morefrom2 ]; then
 echo "merging mli"
 gdal_merge.py -n 0 -co COMPRESS=DEFLATE -o $frame/GEOC/$frame.geo.mli.tif $morefrom1 $morefrom2 
 mv $frame.geo.mli.tif $frame/GEOC/.
fi
fi

for ext in geo.iono.code.tif  sltd.geo.tif  tide.geo.tif geo.mli.tif; do
  echo "merging "$ext
for epoch in `ls $LiCSAR_public/$track/$frame2/epochs | grep 20`; do 
if [ ! -f $frame/epochs/$epoch/$epoch.$ext ]; then
 morefrom1=$LiCSAR_public/$track/$frame2/epochs/$epoch/$epoch.$ext
 morefrom2=$LiCSAR_public/$track/$frame1/epochs/$epoch/$epoch.$ext
 if [ -f $morefrom1 ] && [ -f $morefrom2 ]; then
   mkdir -p $frame/epochs/$epoch
   echo "  "$epoch
   gdal_merge.py -n 0 -co COMPRESS=DEFLATE -o $frame/epochs/$epoch/$epoch.$ext $morefrom1 $morefrom2 # output file name containing path got error
   #mv $epoch.sltd.geo.tif $frame/epochs/$epoch/.
 fi
fi
done
done

for ifg in `cat $frame/ifg_merge_unwrap.list`; do
mkdir -p $frame/GEOC/$ifg

cc1=$frame1/GEOC/$ifg/${ifg}.geo.cc.tif
cc2=$frame2/GEOC/$ifg/${ifg}.geo.cc.tif
if [ ! -f $cc1 ]; then
   cc1=$LiCSAR_public/$track/$frame1/interferograms/$ifg/${ifg}.geo.cc.tif
fi
if [ ! -f $cc2 ]; then
   cc2=$LiCSAR_public/$track/$frame2/interferograms/$ifg/${ifg}.geo.cc.tif
fi
if [ ! -f $cc1 -o ! -f $cc2 ]; then
   echo "Coherence files of interferogram $ifg do not exist! Cancelling..."
   rm -r $frame/GEOC/$ifg
   continue
fi

diff1=$frame1/GEOC/$ifg/${ifg}.geo.diff_pha.tif
diff2=$frame2/GEOC/$ifg/${ifg}.geo.diff_pha.tif
undiff1=$frame1/GEOC/$ifg/${ifg}.geo.diff_unfiltered_pha.tif
undiff2=$frame2/GEOC/$ifg/${ifg}.geo.diff_unfiltered_pha.tif
if [ ! -f $diff1 ]; then
   diff1=$LiCSAR_public/$track/$frame1/interferograms/$ifg/${ifg}.geo.diff_pha.tif
fi
if [ ! -f $diff2 ]; then
   diff2=$LiCSAR_public/$track/$frame2/interferograms/$ifg/${ifg}.geo.diff_pha.tif
fi
if [ ! -f $undiff1 ]; then
   undiff1=$LiCSAR_public/$track/$frame1/interferograms/$ifg/${ifg}.geo.diff_unfiltered_pha.tif
fi
if [ ! -f $undiff2 ]; then
   undiff2=$LiCSAR_public/$track/$frame2/interferograms/$ifg/${ifg}.geo.diff_unfiltered_pha.tif
fi

cur=`cat $frame/ifg_merge_unwrap.list | grep -n $ifg | awk -F ":" '{print $1}'`
echo "Working on interferogram pair $ifg ($cur of $sum) ..."

if [ ! -f $frame/GEOC/$ifg/${ifg}.geo.cc.tif ]; then
   echo "merging cc..."
   gdal_merge.py -n 0 -co COMPRESS=DEFLATE -o $frame/GEOC/$ifg/${ifg}.geo.cc.tif $cc1 $cc2
fi

if [ -f $frame/GEOC/$ifg/${ifg}.geo.unw.tif ]; then
  echo "Interferogram pair $ifg has been merged and unwrapped successfully! Cancelling..."
  continue
else
  unw1=$frame1/GEOC/$ifg/${ifg}.geo.unw.ovlp.tif
  unw2=$frame2/GEOC/$ifg/${ifg}.geo.unw.ovlp.tif

  if [ ! -f $unw1 ]; then
     R=`gmt grdselect $diff1 $diff2 -Ai`
     gmt grdcut $diff1 -G$frame1/GEOC/$ifg/${ifg}.geo.diff_pha.ovlp.tif $R
     gmt grdcut $mask1 -G$frame1/GEOC/${frame1}.geo.landmask.ovlp.tif $R
     gmt grdcut $cc1 -G$frame1/GEOC/$ifg/${ifg}.geo.cc.ovlp.tif $R
     M1=$frame1/GEOC/$ifg/${ifg}.geo.diff_pha.ovlp.tif
     M2=$frame1/GEOC/$ifg/${ifg}.geo.cc.ovlp.tif
     S1=`gmt grdinfo $M1 | grep n_columns | awk '{print $11}'`
     S2=`gmt grdinfo $M2 | grep n_columns | awk '{print $11}'`
     if [ ! $S1 == $S2 ]; then
        gdalwarp2match.py $M2 $M1 $frame1/GEOC/$ifg/${ifg}.geo.cc.ovlp.ok.tif
        M2=$frame1/GEOC/$ifg/${ifg}.geo.cc.ovlp.ok.tif
        gdalwarp2match.py $M1 $M2 $frame1/GEOC/$ifg/${ifg}.geo.diff_pha.ovlp.ok.tif
        M1=$frame1/GEOC/$ifg/${ifg}.geo.diff_pha.ovlp.ok.tif
        mv $M1 $frame1/GEOC/$ifg/${ifg}.geo.diff_pha.ovlp.tif
        mv $M2 $frame1/GEOC/$ifg/${ifg}.geo.cc.ovlp.tif
     fi

     unwrap_geo_ovlp.sh $frame1 $ifg
  fi

  if [ ! -f $unw2 ]; then
     R=`gmt grdselect $diff1 $diff2 -Ai`
     gmt grdcut $diff2 -G$frame2/GEOC/$ifg/${ifg}.geo.diff_pha.ovlp.tif $R
     gmt grdcut $mask2 -G$frame2/GEOC/${frame2}.geo.landmask.ovlp.tif $R
     gmt grdcut $cc2 -G$frame2/GEOC/$ifg/${ifg}.geo.cc.ovlp.tif $R
     M1=$frame2/GEOC/$ifg/${ifg}.geo.diff_pha.ovlp.tif
     M2=$frame2/GEOC/$ifg/${ifg}.geo.cc.ovlp.tif
     S1=`gmt grdinfo $M1 | grep n_columns | awk '{print $11}'`
     S2=`gmt grdinfo $M2 | grep n_columns | awk '{print $11}'`
     if [ ! $S1 == $S2 ]; then
        gdalwarp2match.py $M2 $M1 $frame2/GEOC/$ifg/${ifg}.geo.cc.ovlp.ok.tif
        M2=$frame2/GEOC/$ifg/${ifg}.geo.cc.ovlp.ok.tif
        gdalwarp2match.py $M1 $M2 $frame2/GEOC/$ifg/${ifg}.geo.diff_pha.ovlp.ok.tif
        M1=$frame2/GEOC/$ifg/${ifg}.geo.diff_pha.ovlp.ok.tif
        mv $M1 $frame2/GEOC/$ifg/${ifg}.geo.diff_pha.ovlp.tif
        mv $M2 $frame2/GEOC/$ifg/${ifg}.geo.cc.ovlp.tif
     fi

     unwrap_geo_ovlp.sh $frame2 $ifg
  fi

  if [ -f $unw1 -a -f $unw2 ]; then
     gdalwarp2match.py $unw2 $unw1 $frame2/GEOC/$ifg/${ifg}.unw.ok.tif
     unw2=$frame2/GEOC/$ifg/${ifg}.unw.ok.tif
     gdalwarp2match.py $unw1 $unw2 $frame1/GEOC/$ifg/${ifg}.unw.ok.tif
     unw1=$frame1/GEOC/$ifg/${ifg}.unw.ok.tif
     # seems grdmath must always save to some grd, pity..
     gmt grdmath $unw2 0 NAN $unw1 0 NAN SUB MEDIAN = $frame/GEOC/$ifg/median.grd
     correction=`gmt grdinfo -T $frame/GEOC/$ifg/median.grd | cut -d '/' -f2`
     gmt grdmath $diff1 0 NAN $correction ADD WRAP = $frame1/GEOC/$ifg/${ifg}.geo.diff_pha.corr.grd
     gmt grdmath $undiff1 0 NAN $correction ADD WRAP = $frame1/GEOC/$ifg/${ifg}.geo.diff_unfiltered_pha.corr.grd
     gdal_translate -of GTiff -ot Float32 -co COMPRESS=DEFLATE -a_srs EPSG:4326 $frame1/GEOC/$ifg/${ifg}.geo.diff_pha.corr.grd $frame1/GEOC/$ifg/${ifg}.geo.diff_pha.corr.tif
     gdal_translate -of GTiff -ot Float32 -co COMPRESS=DEFLATE -a_srs EPSG:4326 $frame1/GEOC/$ifg/${ifg}.geo.diff_unfiltered_pha.corr.grd $frame1/GEOC/$ifg/${ifg}.geo.diff_unfiltered_pha.corr.tif
     gdal_merge.py -n 0 -co PREDICTOR=3 -co COMPRESS=DEFLATE -o $frame/GEOC/$ifg/${ifg}.geo.diff_pha.tif $frame1/GEOC/$ifg/${ifg}.geo.diff_pha.corr.tif $diff2
     gdal_merge.py -n 0 -co PREDICTOR=3 -co COMPRESS=DEFLATE -o $frame/GEOC/$ifg/${ifg}.geo.diff_unfiltered_pha.tif $frame1/GEOC/$ifg/${ifg}.geo.diff_unfiltered_pha.corr.tif $undiff2
     # diff_unfiltered_pha
     #mv ${ifg}.geo.diff_pha.tif $frame/GEOC/$ifg/
 
     if [ $LOTUS -eq 1 ]; then 
      bsub2slurm.sh -o $frame'_'$ifg.out -e $frame'_'$ifg.err -J 'mergeunw_'$frame'_'$ifg -q comet -n 1 -W 00:50 -M 12288 jin_unwrap_geo.sh $frame $ifg
     else
      jin_unwrap_geo.sh $frame $ifg
     fi
  else
     echo $ifg >> $frame/ifg_failed.list
  fi
fi
done

if [ -f $frame/ifg_failed.list ]; then
   sum=`cat $frame/ifg_failed.list | wc -l`
   echo "$sum interferogram pairs failed... further details in $frame/ifg_failed.list"
   #cat $frame/ifg_failed.list
else
   echo "$sum interferogram pairs have been successfully merged and unwrapped!"
fi

echo "correcting directory structure"
mkdir $frame/metadata
mv $frame/GEOC/*tif $frame/metadata/.
cp $frame1/metadata/metadata.txt $frame/metadata/.
mv $frame/GEOC $frame/interferograms
stop_time=`date +%s`
echo "Elapsed time: `echo "scale=1;($stop_time - $start_time)/3600" | bc` h"

