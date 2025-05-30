#! /bin/bash
module load ICAMS
source $LiCSARpath/lib/LiCSAR_bash_lib.sh  

if [ -z $1 ]; then
 echo "Usage: icams_one_epoch.sh FRAMEID EPOCHDATE [CLEANIT]"
 exit
fi


frame=$1
epoch=$2
cleanit=0
if [ ! -z $3 ]; then
  cleanit=$3
fi

framedir=$LiCSAR_public/`track_from_frame $frame`/$frame
epochsdir=$framedir/epochs
metadir=$framedir/metadata
epochdir=$epochsdir/$epoch

icamsout=icams/ERA5/sar/$epoch'_tot.h5'
sltdout=$epochdir/$epoch.icams.sltd.geo.tif

if [ -f $sltdout ]; then
  echo "ICAMS correction already exists, cancelling for epoch "$epoch
  exit
fi

if [ ! -d $epochdir ]; then mkdir $epochdir; fi

hgt=$metadir/$frame.geo.hgt.tif
U=$metadir/$frame.geo.U.tif

if [ ! -f $U ]; then echo "ERROR, U file not exists, cancelling"; exit; fi

wesn=`gmt grdinfo -Ir $hgt | cut -d 'R' -f2`
resdeg=`gmt grdinfo -I $hgt | cut -d '/' -f2`
resol=`echo $resdeg" * 111111" | bc | cut -d '.' -f1`

ttime=`grep center_time $metadir/metadata.txt | cut -d '=' -f2 | cut -c -5`

#ppwd=`pwd`
#epoch=`basename $ppwd`
### FIX for frames thata SLC is more 23:30 P. Espin

echo $ttime
# Convert ttime to minutes since midnight for comparison
ttime_minutes=$(echo $ttime | awk -F: '{ print ($1 * 60) + $2 }')
limit_minutes=$((23 * 60 + 30))
# Combine date and time into a single variable
datetime="${date} ${ttime}"

if (( ttime_minutes > limit_minutes )); then
    tt="00:00"
    # Increment the date by one day and format it back to YYYYMMDD
    # Convert input date to a recognized format for date manipulation
    formatted_date=$(date -d "$epoch" +%Y-%m-%d)
    # Increment the date by one day and format it back to YYYYMMDD
    echo $formatted_date
    new_date=$epoch #(date -d "$formatted_date +1 day" +%Y%m%d)    
    echo "New date" ${new_date}
    echo $tt
    icamsout=icams/ERA5/sar/${new_date}'_tot.h5'

else
    tt="$ttime"
    new_date=$epoch
    echo $new_date
    echo $tt
    icamsout=icams/ERA5/sar/$epoch'_tot.h5'

fi

python tropo_icams_date.py ${new_date} --region " "$wesn" " --imaging-time $tt --dem-tif $hgt --resolution $resol

#####


#tropo_icams_date.py $epoch --region " "$wesn" " --imaging-time $ttime --dem-tif $hgt --resolution $resol

python3 -c "import lics_processing as lp; lp.ztd2sltd('"$icamsout"', '"$U"', outif = '"$sltdout"')"

if [ $cleanit -gt 0 ]; then
# clean?
rm -r icams ttt
fi

