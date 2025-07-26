#!/bin/bash

frame=$1
prod_dir=$LiCSAR_public
# old gacos dir:
#work_dir="/gws/nopw/j04/nceo_geohazards_vol2/LiCS/temp/GACOS"
work_dir=$LiCSAR_GACOS
if [ -z $work_dir ]; then
   work_dir='/work/scratch-pw3/licsar/GACOS'
fi

if [ -z $frame ]; then
 echo "Usage: parameter is frame"
 exit
fi


track=`echo $frame | cut -c -3 | sed 's/^0//' | sed 's/^0//'`

if [ -f $prod_dir/$track/$frame/metadata/$frame.geo.U.tif ]; then
   cp $prod_dir/$track/$frame/metadata/$frame.geo.U.tif $work_dir/$frame/.
else
   echo "File  does not exist."
   echo $prod_dir/$track/$frame/metadata/$frame.geo.U.tif
   exit 1
fi

for rsc in $work_dir/$frame/*.ztd; do 
  if [ `ls -al $rsc | gawk {'print $5'}` -lt 100 ]; then rm $rsc; else
   if [ ! -f $rsc.geo.tif ]; then
    LiCSAR_ztd2geotiff.py -z $rsc -u $work_dir/$frame/$frame.geo.U.tif
    # but sometimes it again fails, creating bad (empty) file
    if [ `ls -al $rsc.geo.tif | gawk {'print $5'}` -lt 100 ]; then
      echo "error processing "$rsc
      rm $rsc.geo.tif
    fi
   fi
  fi
done

for geotif in $work_dir/$frame/*.sltd.geo.tif; do
   #date=$(echo $geotif | awk -F"." '{print $1}' | awk -F'/' '{print $12}')
   date=`basename $geotif | cut -d '.' -f1`
   mkdir -p  $prod_dir/$track/$frame/epochs/$date
   cp -f $geotif $prod_dir/$track/$frame/epochs/$date/.
   touch $geotif.copdem
   cp -f $geotif.copdem $prod_dir/$track/$frame/epochs/$date/.
   chmod -R 775 $prod_dir/$track/$frame/epochs/$date 2>/dev/null
   cedaarch_create_html.sh $frame $date
done

#for jpg in $work_dir/$frame/*.ztd.jpg; do
   #date=$(echo $jpg | awk -F"." '{print $1}' | awk -F'/' '{print $12}')
#   date=`basename $jpg | cut -d '.' -f1`
#   cp -f $jpg $prod_dir/$track/$frame/epochs/$date/.
#done


