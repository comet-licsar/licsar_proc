#!/bin/bash
if [ -z $1 ]; then 
 echo "just include ifg pair"
 exit
fi
pair=$1
infile=IFG/$pair/$pair.filt.diff
tmpfile=IFG/$pair/tmpfile
master=`ls geo/20??????.hgt | cut -d '.' -f1 | cut -d '/' -f2 | head -n1`
parfile=SLC/$master/$master.slc.mli.par
octavepath=$LiCSARpath/octave

if [ -f $infile.backup ]; then
 echo "backup file already exists. cancelling as this means we already have run this"
 exit
fi

NOLINES=`grep "azimuth_lines:" $parfile | gawk '{print $2}'`

if [ ! -z `which byteswap 2>/dev/null` ]; then
 byteswap -o $tmpfile 4 $infile
else
 swap_bytes $infile $tmpfile 4
fi

echo "addpath(genpath('"$octavepath"'));" > detrend.m
echo "a=freadbk('"$tmpfile"',"$NOLINES",'cpxfloat32');" >> detrend.m
echo "b=cpxdetrend(a);" >> detrend.m
echo "fwritebk(b,'"$tmpfile"','cpxfloat32');" >> detrend.m
#echo "exit;" >> detrend.m

#if [ ! -z `which matlab 2>/dev/null` ]; then
# cat detrend.m | matlab -nodesktop -nosplash -nojvm -nodisplay
#else
octave-cli -q detrend.m
#fi

mv $infile $infile.backup

if [ ! -z `which byteswap 2>/dev/null` ]; then
  byteswap -o $infile 4 $tmpfile
else
  swap_bytes $tmpfile $infile 4
fi

rm $tmpfile detrend.m 2>/dev/null
