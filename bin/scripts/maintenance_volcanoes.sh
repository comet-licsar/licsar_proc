#!/bin/bash

cd $LiCSAR_procdir/../../volc-proc/current
for area in `ls`; do
  echo $area
  for fld in `ls $area`; do
    echo $fld
    fr=`echo $fld | rev | cut -c -17 | rev`
    isz=0
    if [ -d $area/$fld/tif ]; then
      for fl in `ls $area/$fld/tif`; do
        if [ ! -s $area/$fld/tif/$fl ]; then
          isz=1
          rm $area/$fld/tif/$fl
        fi
      done
      if [ $isz -gt 0 ]; then
        falbino_volc_clip_figs.py $area $fr
      fi
    fi
  done
done