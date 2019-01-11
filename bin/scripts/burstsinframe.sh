#!/bin/bash 
#################################################################
# This program identify which bursts are in 
# a given directory containing multiple zipfiles
#
#  Dependencies:
#   extract_bursts_loc_v3.0.py 
#   BurstList_2_slcID.py
#
# Author: Richard Walters 
# Date: 2016/03/10
# Modified: add documentation and re-code for stable version (PJG)
#################################################################

usage() { 
  echo " " 1>&2; 
  echo "Usage: $0 " 1>&2; 
  echo "   This program identify which bursts are in a given directory containing multiple zipfiles  " 1>&2; 
  echo "   e.g. burstsinframe.sh A14frame3.xy A14frame3_burst_ids.txt /nfs/a1/raw/sentinel/T014A/20151130 151130 14 " 1>&2; 
  echo " " 1>&2; 
  exit 1; 
}

if [ "$#" -gt 5 ]; then
  usage
else
  polyfile=$1
  burstids=$2
  SLCdir=$3
  date=$4
  track=$5
fi

export PATH=/nfs/a1/insar/sentinel1/lics_bin/Burst_stuff/:$PATH
rm -f ${date}_zipfiles.txt

# extract bursts from each SLC, cat together
# remove previous and create new blank files
rm -f allcenterIW1.list allcenterIW2.list allcenterIW3.list
touch allcenterIW1.list allcenterIW2.list allcenterIW3.list 

ls -d $SLCdir/S1*$date*.zip > tmp_ziplist.txt
awk -F"/" '{a  = substr($NF,50,6)%175 -72; if (a < 0) track = a+175; else track = a; if (track == '$track') print $0}' tmp_ziplist.txt > tmp_ziplist2.txt
for file in `more tmp_ziplist2.txt`; do
    extract_bursts_loc_v3.0.py $file
    #cat together coordcenters
    awk '{print $1, $2, $3, 1, "'$file'"}' coordcenterIW1.list > coord1names.list
    cat allcenterIW1.list coord1names.list > tmp1.txt
    mv tmp1.txt allcenterIW1.list
    awk '{print $1, $2, $3, 2, "'$file'"}' coordcenterIW2.list > coord2names.list
    cat allcenterIW2.list coord2names.list > tmp2.txt
    mv tmp2.txt allcenterIW2.list
    awk '{print $1, $2, $3, 3, "'$file'"}' coordcenterIW3.list > coord3names.list
    cat allcenterIW3.list coord3names.list > tmp3.txt
    mv tmp3.txt allcenterIW3.list
    rm coord1names.list coord2names.list coord3names.list
done

cat allcenterIW1.list allcenterIW2.list allcenterIW3.list > allcenter.list
awk '{print $1, $2}' allcenter.list > lonlats.tmp
listlonlat2mgrs.py lonlats.tmp > tmpids
paste tmpids allcenter.list > tmp.txt 
mv tmp.txt allcenter.list

rm -f tmp1 tmp2 ${date}_zipfiles.txt

nbursts=`wc -l $burstids | awk '{print $1 -1}'`
for i in `seq 1 $nbursts`; do
    id=`awk '(NR==('$i' +1)){print $1}' $burstids`
    lon=`awk '(NR==('$i' +1)){print $2}' $burstids`
    lat=`awk '(NR==('$i' +1)){print $3}' $burstids`
    awk '{print $0, 111*sqrt((($2-'$lon')*cos('$lat'*(3.14159/180)))^2 + ($3-'$lat')^2)}' allcenter.list > tmplist
    match=`sort -k7 -n tmplist | head -1`
    echo $id $lon $lat $match >> ${SLCdir}/${date}_zipfiles.txt
done

awk '{if ($10<5) print $2, $3, $7, $8, $9, $10; else print  $2, $3, $7, $8, "NA", $10}' ${SLCdir}/${date}_zipfiles.txt > ${SLCdir}/${date}_zipfiles_short.txt


rm -f allcenter*.list coord*.list lonlats.tmp tmp* 
