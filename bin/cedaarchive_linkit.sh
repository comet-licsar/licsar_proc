#!/bin/bash

LiCSAR_public=/gws/nopw/j04/nceo_geohazards_vol1/public/LiCSAR_products_test
cd $LiCSAR_public
frame=041A_04785_131313
track=`track_from_frame $frame`
mkdir -p $track/$frame/interferograms
cd $track/$frame/interferograms

for x in `ls /neodc/comet/data/licsar_products/$track/$frame/*_* -d`; do 
 pair=`basename $x`
 for fname in `ls $x`; do
  echo "<a href='https://data.ceda.ac.uk/neodc/comet/data/licsar_products/"$track"/"$frame"/"$pair"/"$fname"'>"$fname"</a><br />" >> $pair;
 done
done
