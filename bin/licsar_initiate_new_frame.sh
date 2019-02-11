#!/bin/bash
#a script that will generate master, DEM  and all other necessary files to initiate a new frame
#module load licsar_proc
curdir=$LiCSAR_procdir

if [ -z $1 ];
then
 echo "Usage: licsar_initiate_new_frame.sh $FRAME"
 echo "where frame can be e.g. 010D_11058_131313"
 exit
else
 frame=$1
 tr=`echo $frame | cut -d '_' -f1 | sed 's/^0//' | sed 's/^0//' | rev | cut -c 2- | rev`
 if [ -d $curdir/$tr/$frame ]; then
  echo "This frame already exists! Stopping here"
  echo "Check and remove(?) "$curdir/$tr/$frame
  exit
 fi
 if [ `echo $frame | grep -o '_' | wc -l` != 2 ]; then
  echo "Wrong frame name. Stopping"
  exit
 fi
fi

echo "Setting the master image and DEM for frame "$frame
LiCSAR_setup_master.py -f $frame -d $curdir/$tr/$frame -A -r 20 -a 4
cd $curdir/$tr/$frame
LiCSAR_05_mk_angles_master

echo "cleaning"
rm -f $curdir/tr/$frame/SLC/*/2*T*.I*sl* 2>/dev/null

echo "done"
