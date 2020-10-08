#!/bin/bash

if [ -z $1 ]; then echo "parameter is FRAME"; exit; fi

frame=$1
track=`echo $frame | cut -c -3 | sed 's/^0//' | sed 's/^0//'`
outfile=$LiCSAR_public/$track/$frame/metadata/metadata.txt

master=`ls $LiCSAR_procdir/$track/$frame/geo/*[0-9].hgt | rev | cut -d '/' -f1 | rev | cut -d '.' -f1 | head -n1`

#write master time to the outfile
echo "master="$master > $outfile

#write coarse acquisition time (e.g. for GACOS purposes)
time=`grep center_time $LiCSAR_procdir/$track/$frame/SLC/$master/$master.slc.par | gawk {'print $2'}`
time_hours=`python -c "import datetime as dt; time=dt.timedelta(seconds="$time"); print(str(time))"`
echo "center_time="$time_hours >> $outfile
