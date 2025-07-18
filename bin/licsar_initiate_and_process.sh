#!/bin/bash
#curdir=$LiCSAR_procdir

if [ -z $1 ];
then
 echo "Usage: licsar_initiate_and_process.sh FRAME REFDATE 2016-02-01 2024-06-01"
 echo "this will send frame init to LOTUS2, and run background update for the dates range"
 echo "to know REFDATE (in the form of yyyymmdd), you can just run it with only FRAME.."
 exit
fi

frame=$1

if [ -z $2 ]; then
  licsar_initiate_new_frame.sh -V 2000 $frame
  echo "please pick up some and set as m=yyyymmdd"
  exit
fi

m=$2
startdate=$3
enddate=$4

#tr=`track_from_frame $fr`
#relorb=`echo $fr | cut -d '_' -f 1`
#get_dates_scihub.py $fr 20180101 20210101 | grep 20 > $fr.mdates
#ls $LiCSAR_procdir/$tr/$relorb'_'*/SLC | grep ^20 > $fr.surrounding
#echo "best options perhaps:"
#for x in cat $fr.surrounding; do grep $x $fr.mdates >> $fr.mcands; done
#cat $fr.mcands
#for m in $fr.mcands; do
#for m in cat $fr.mdates; do
# echo $m
# check_master_candidate.py $fr $m
#done

#m=20190729
echo licsar_initiate_new_frame.sh $fr $m > $fr.init
chmod 777 $fr.init
bsub2slurm.sh -o $fr.init.out -e $fr.init.err -J $fr.init -n 1 -W 10:59 -M 32678 ./$fr.init
nohup framebatch_update_frame.sh $fr gapfill 2016-02-01 2024-06-01 >$fr.init.proc.out 2>$fr.init.proc.err &

