#!/bin/bash

frame=$1
tr=`track_from_frame $frame`

if [ ! -d $LiCSAR_procdir/$tr/$frame ]; then
  echo "error - this frame has no directory in LiCSAR_procdir..."
  exit
fi

for s in `ls $LiCSAR_procdir/$tr/$frame/subsets -alh 2>/dev/null | grep subsets | gawk {'print $11'}`; do
  mv $s $s.backup
done

for x in $LiCSAR_public $LiCSAR_procdir $LiCSAR_web; do
  if [ -d $x/$tr/$frame.backup ]; then echo "error - this backup exists"; exit; fi
done

for x in $LiCSAR_public $LiCSAR_procdir $LiCSAR_web; do
  mv $x/$tr/$frame $x/$tr/$frame.backup
done

mv $BATCH_CACHE_DIR/$frame $BATCH_CACHE_DIR/$frame.backup 2>/dev/null

echo "the frame was moved to "$frame".backup directories"