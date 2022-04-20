#!/bin/bash

cd $LiCSAR_public
for tr in `seq 1 175`; do
echo $tr
for fr in `ls $tr`; do
for ifg in `ls $tr/$fr/interferograms 2>/dev/null`; do

for meta in `ls $tr/$fr/interferograms/$ifg/*metadata 2>/dev/null`; do
rm $meta
todo=`basename $meta .metadata`
epos_export_tif.sh $LiCSAR_public/$todo
done

for kmz in `ls $tr/$fr/interferograms/$ifg/*kmz 2>/dev/null`; do
for todoext in cc diff_pha unw; do
todo=$tr/$fr/interferograms/$ifg/$ifg.geo.$todoext.tif
if [ ! -f $todo.metadata ]; then
epos_export_tif.sh $LiCSAR_public/$todo
fi
done
done

done
done
done

