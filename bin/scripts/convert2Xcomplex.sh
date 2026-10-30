#!/bin/bash

if [ -z $2 ]; then
  echo "Usage: path to SLC or RSLC epoch dir - this will change it to either SCOMPLEX or FCOMPLEX"
  echo "e.g. convert2Xcopmlex.sh ./RSLC/20160914 FCOMPLEX"
  echo "(it will first check if it is needed here, and will fix also if it is in both SLC and RSLC - e.g. ref epoch)"
  exit
fi

ctype=$2
if [ ! `echo $ctype | cut -c 2-` == 'COMPLEX' ]; then
  echo "Wrong type - use either SCOMPLEX or FCOMPLEX"
  exit
fi
if [ $ctype == 'FCOMPLEX' ]; then origctype='SCOMPLEX'; else origctype='FCOMPLEX'; fi

epochpath=$1
epoch=`basename $epochpath`
for slc in `ls $epochpath/../../*LC/$epoch/*slc`; do
  slcpar=$slc.par
  if [ -f $slcpar ]; then
    # first check if it is link...
    if [ -L $slc ]; then
      # maybe it was already changed?
      realslc=`realpath $slc`
      realslcpar=$realslc.par
      if [ `grep -c $ctype $realslcpar` -lt 1 ]; then
        # if not, it means, the real datafile was already changed ..
        echo "converting "$slc
        rslc2rslc.py $slc $slc.xcx 0 $ctype; rm $slc; mv $slc.xcx $slc
      fi
    elif [ `grep -c $ctype $slcpar` -lt 1 ]; then
      echo "converting "$slc
      rslc2rslc.py $slc $slc.xcx 0 $ctype; rm $slc; mv $slc.xcx $slc
    fi
    sed -i 's/'$origctype'/'$ctype'/' $slcpar
  fi
done