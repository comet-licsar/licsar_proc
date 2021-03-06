#!/bin/bash
#a script that will generate master, DEM  and all other necessary files to initiate a new frame
#module load licsar_proc #or _testing
if [ `basename $LiCSAR_procdir` != "eq" ]; then LiCSAR_procdir=$LiCSAR_procdir/eq; fi
curdir=$LiCSAR_procdir

if [ -z $1 ];
then
 echo "Usage: licsar_initiate_new_frame_eq.sh FRAME [MASTER_YYMMDD]"
 echo "where frame can be e.g. 010D_11058_131313 20180722"
 echo "(if master date is not given, it will choose automatically, from last 3 months data)"
 echo "(to include custom downloaded master files, don't forget to arch2DB.py them first)"
 exit
fi
if [ ! -z $2 ]; then
 getmaster="-m "$2
 else
 getmaster="-A"
fi

 frame=$1
 tr=`echo $frame | cut -d '_' -f1 | sed 's/^0//' | sed 's/^0//' | rev | cut -c 2- | rev`
 mkdir $curdir/$tr 2>/dev/null
 #if [ -d $curdir/$tr/$frame ]; then
 # echo "This frame already exists! Stopping here"
 # echo "Check and remove(?) "$curdir/$tr/$frame
 # exit
 #fi
 if [ `echo $frame | grep -o '_' | wc -l` != 2 ]; then
  echo "Wrong frame name. Stopping"
  exit
 fi

echo "Setting the master image and DEM for frame "$frame
LiCSAR_setup_master.py -f $frame -d $curdir/$tr/$frame $getmaster -r 20 -a 4 -e 1
if [ ! -d $curdir/$tr/$frame/SLC ]; then
 echo "Something got wrong with the initiation"
else
 cd $curdir/$tr/$frame
 m=`ls SLC`
 if [ ! -d RSLC/$m ]; then 
  mkdir RSLC/$m
  for slcfile in `ls $curdir/$tr/$frame/SLC/$m/*`; do ln -s $slcfile `echo $slcfile | sed 's/SLC/RSLC/' | sed 's/slc/rslc/'`; done
 fi
# LiCSAR_05_mk_angles_master
 #echo "Generating E-N-U files"
 #submit_lookangles.py -f $frame -t $tr
 #sometimes .xy is not generated..
 if [ ! -f $curdir/$tr/$frame/frame.xy ]; then cp $curdir/$tr/$frame/$frame'-poly.txt' $curdir/$tr/$frame/frame.xy; fi
 #cp $curdir/$tr/$frame/$frame'-poly.txt' $LiCSAR_public/$tr/$frame/metadata/.
 echo "cleaning"
 rm -f $curdir/$tr/$frame/SLC/*/2*T*.I*sl* 2>/dev/null
 echo "done"
fi
