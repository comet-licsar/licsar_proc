#!/bin/bash

if [ -z $2 ]; then
 echo "usage: remove_from_lics.sh FRAMEID EPOCH (or IFG) OPTIONAL"
 echo "e.g. remove_from_lics.sh 144D_05298_131313 20190907 0"
 echo "OPTIONAL=0 delete the subset clip for high resolution volcanoes RSLC"
 exit
fi

if [ "$#" -lt 3 ]; then
    echo "Uso: $0 <param1> <param2> <param3> [<opcional>]"
    exit 1
fi


frame=$1
epoch=$2
optional=$3


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

master=`ls $procdir/SLC | head -n1`

if [ $epoch == $master ]; then
 echo "this is a master file - cannot delete it"
 exit
fi
#echo $pubdir
#echo $procdir
#echo 'bbbb'$epoch'bbbb'
#echo "ls -d $pubdir/*/*$epoch*"

for dir in $pubdir/epochs $pubdir/interferograms $procdir/log $procdir/tab $procdir/LUT $procdir/RSLC; do
 echo "deleting:"
 ls -d $dir/*$epoch* 2>/dev/null
 rm -rf $dir/*$epoch* 2>/dev/null
done



# Si se incluye el cuarto argumento y es "0", imprime "SI"
if [ "$optional" == "0" ]; then
    echo "SI"


    curdir=$LiCSAR_procdir
    tr=`echo $frame | cut -d '_' -f1 | sed 's/^0//' | sed 's/^0//' | rev | cut -c 2- | rev`
    tr1=`(echo $frame | cut -d '_' -f1 | sed 's/^0*//' | sed 's/[A-Z]$/D/')`
    echo $tr1
    frameDir=$curdir/$tr/$frame


    if [ -d $frameDir/subsets ]; then
       echo "clipping for subsets"
       for subset in `ls $frameDir/subsets`; do
         echo "subset "$subset
         cornersclip=$frameDir/subsets/$subset/corners_clip.$frame
         subdir=$frameDir/subsets/$subset
         dirsub="/gws/nopw/j04/nceo_geohazards_vol1/projects/LiCS/proc/current/$subset/volc/$tr1/RSLC/$epoch"
         #echo $dirsub
         #echo $subset
         echo "delete $epoch from $subset"

       done
    fi
fi


