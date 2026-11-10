#!/bin/bash

if [ -z $2 ]; then
  echo "Usage: path to SLC or RSLC epoch dir - this will change it to either SCOMPLEX or FCOMPLEX"
  echo "e.g. convert2Xcomplex.sh ./RSLC/20160914 FCOMPLEX"
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
ftype=`ls $epochpath/$epoch.*slc | head -n 1 | rev | cut -d '.' -f 1 | rev`  # rslc or slc
if [ -z $ftype ]; then echo "ERROR - seems no (r)slc file exists in the directory"; exit; fi
doit=0
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
        doit=1
        #echo "converting "$slc
        #rslc2rslc.py $slc $slc.xcx 0 $ctype; rm $slc; mv $slc.xcx $slc
      fi
    elif [ `grep -c $ctype $slcpar` -lt 1 ]; then
      #echo "converting "$slc
      doit=1
      #rslc2rslc.py $slc $slc.xcx 0 $ctype; rm $slc; mv $slc.xcx $slc
    fi
    # sed -i 's/'$origctype'/'$ctype'/' $slcpar
  fi
done

# updated version to use GAMMA script
if [ $doit -gt 0 ]; then
    echo "converting "$epoch
    mkdir -p SLC.temp.fcpx tab
    mv $epochpath SLC.temp.fcpx/.; mkdir $epochpath
    createSLCtab SLC.temp.fcpx/$epoch/$epoch $ftype > tab/$epoch.fcpx.tab
    createSLCtab $epochpath/$epoch $ftype > tab/$epoch.tab
    rm -f tab/$epoch.burstlist
    for x in SLC.temp.fcpx/$epoch/$epoch.IW?.$ftype.TOPS_par; do
      echo "1 "`get_value $x number_of_bursts` >> tab/$epoch.burstlist
    done
    echo "converting "$epoch" to SCOMPLEX"
    SLC_copy_ScanSAR tab/$epoch.fcpx.tab tab/$epoch.tab tab/$epoch.burstlist 1 >/dev/null 2>/dev/null
    if [ $ftype == 'slc' ]; then
      echo "mosaicing (maybe not needed)"
      SLC_mosaic_ScanSAR tab/$epoch.tab $epochpath/$epoch.slc $epochpath/$epoch.slc.par 20 4 1 >/dev/null 2>/dev/null
      multilookSLC $epoch 20 4
    fi
    rm -r SLC.temp.fcpx/$epoch tab/$epoch.fcpx.tab
    rmdir SLC.temp.fcpx 2>/dev/null
else
  echo "nothing to do"
fi