#!/bin/bash

if [ -z $2 ]; then
 echo "usage: remove_from_lics.sh FRAMEID EPOCH (or IFG)"
 echo "e.g. remove_from_lics.sh 144D_05298_131313 20190907"
 exit
fi

frame=$1
epoch=$2

track=`echo $frame | cut -c -3 | sed 's/^0//' | sed 's/^0//'`
#epoch=`echo $epoch | sed 's/ //g'`

if [ `echo $epoch | grep -c '_'` -eq 1 ]; then
 #ifg
 epoch=`echo $epoch | cut -c -17`
else
 #epoch
 epoch=`echo $epoch | cut -c -8`
fi
if [ ! `echo $epoch | cut -c -2` -eq 20 ]; then
 echo "wrong second parameter (epoch or ifg to delete)"
 exit
fi

changefile=/gws/nopw/j04/nceo_geohazards_vol1/public/LiCSAR_products/frameid_changes.txt
if [ `grep -c ^$frame $changefile` -gt 0 ]; then
 line=`grep ^$frame $changefile`;
 echo "frame ID changed: " $line;
 old=`echo $line | cut -d ',' -f1`;
 new=`echo $line | cut -d ',' -f2`;
 frame=$new
fi

pubdir=$LiCSAR_public/$track/$frame
procdir=$LiCSAR_procdir/$track/$frame
batchfrdir=$BATCH_CACHE_DIR/$frame
webdir=$LiCSAR_web/$track/$frame

master=`ls $procdir/SLC | head -n1`

if [ $epoch == $master ]; then
 echo "this is a master epoch - cannot delete it"
 exit
fi


# if this is epoch, it is probably bad one, so will remove its SD values:
if [ `echo $epoch | wc -m` == 9 ]; then
 python3 -c "import LiCSquery as lq; lq.delete_esds_for_frame('"$frame"', epoch = '"$epoch"', test=False)"
 # also for epoch only - remove from subsets
 if [ -d $procdir/subsets ]; then
  echo "deleting from subsets (if epoch)"
  rm -rf $procdir/subsets/*/RSLC/$epoch #2>/dev/null
 fi
fi

#echo $pubdir
#echo $procdir
#echo 'bbbb'$epoch'bbbb'
#echo "ls -d $pubdir/*/*$epoch*"

# adding to log file
logfile=/gws/nopw/j04/nceo_geohazards_vol1/public/LiCSAR_products/updates/`date +'%Y%m%d'`.$frame.removed
touch $logfile; chmod 775 $logfile
find $pubdir -name '*'$epoch'*' >> $logfile

for dir in $pubdir $procdir $batchfrdir $webdir; do
 echo "deleting:"
 ls -d $dir/*/*$epoch* 2>/dev/null
 rm -rf $dir/*/*$epoch*
done

