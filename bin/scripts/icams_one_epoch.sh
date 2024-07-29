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

tropo_icams_date.py $epoch --region " "$wesn" " --imaging-time $ttime --dem-tif $hgt --resolution $resol

python3 -c "import lics_processing as lp; lp.ztd2sltd('"$icamsout"', '"$U"', outif = '"$sltdout"')"

if [ $cleanit -gt 0 ]; then
# clean?
rm -r icams ttt
fi

