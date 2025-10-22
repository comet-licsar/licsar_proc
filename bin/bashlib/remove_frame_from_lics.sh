#!/bin/bash

frame=$1
tr=`track_from_frame $frame`

for x in $LiCSAR_public $LiCSAR_procdir $LiCSAR_web; do
  if [ -d $x/$tr/$frame.backup ]; then echo "error - this backup exists"; exit; fi
done

for x in $LiCSAR_public $LiCSAR_procdir $LiCSAR_web; do
  mv $x/$tr/$frame $x/$tr/$frame.backup
done
