#!/bin/bash

#################################################################
# This program queries sci hub to find zip files which are within a polygon.
# Though optional, one should supply a track number.
#
#  Dependencies:
#
# Author: David Mackenzie
# Date: 2016/11/04
#################################################################

usage() { 
  echo " " 1>&2; 
  echo "  This program queries sci hub to find zip files which are within a polygon. " 1>&2; 
  echo "Usage: " 1>&2; 
  echo "  query_sentinel.sh track_no asc/dsc polyfile.xy startdate stopdate " 1>&2; 
  echo "   " 1>&2; 
  echo "  " 1>&2; 
  echo "  - can be used instead to collect all tracks within polygon (BETA feature!)" 1>&2; 
  echo "  " 1>&2; 
  echo "  Either asc => ascending orbit or dsc => descending orbit" 1>&2; 
  echo "  " 1>&2; 
  echo "  - can be used instead to collect both ascending and descending (BETA feature!)" 1>&2; 
  echo "  " 1>&2; 
  echo "  start/stop date: yyyymmdd" 1>&2; 
  echo "  " 1>&2; 
#  echo "  Authentication: Is done using your scihub login. Set these details in your .netrc file" 1>&2; 
#  echo "  		e.g. touch ~/.netrc; chmod 600 ~/.netrc; " 1>&2; 
#  echo "  		echo machine scihub.copernicus.eu login MyUserName password MyPassword >> ~/.netrc" 1>&2; 
#  echo "  " 1>&2; 
  echo "  Output: Two files will be created in the same location as the polygon file" 1>&2; 
  echo "          - polyfile_zipfiles.list " 1>&2; 
  echo "          - polyfile_urls.list " 1>&2; 
  echo "  " 1>&2; 
  exit 1; 
}

if [ "$#" -ne 5 ]; then
  usage
else
  track=$1
  ascdsc=$2
  polyfile=$3
  startdate=$4
  stopdate=$5
fi

#getting licsar username/password for scihub
export WGETRC=$LiCSAR_configpath/scihub_credentials_wget

stub=`echo $polyfile | sed 's/\.xy//'`
ziplist=${stub}_zipfile_names.list
linklist=${stub}_urls.list

echo "Track number: $track"
echo "Orbit direction: $ascdsc"
echo "Polygon: $polyfile"
echo "Start epoch: $startdate"
echo "End epoch: $stopdate"

nmax=100

# Parse ascending/descending arguemnt.
if [ $ascdsc == "asc" ]; then
  orbitdir="ASCENDING"
elif [  $ascdsc == "dsc" ]; then
  orbitdir="DESCENDING"
elif [  $ascdsc == "-" ]; then
  orbitdir='-'
else
  echo "Bad orbit direction. Enter asc or dsc."
  exit 1;
fi


# Split dates into components:
startyyyy=`echo $startdate | awk 'BEGIN{FIELDWIDTHS="4 2 2"}{print $1}'`
startmm=`echo $startdate | awk 'BEGIN{FIELDWIDTHS="4 2 2"}{print $2}'`
startdd=`echo $startdate | awk 'BEGIN{FIELDWIDTHS="4 2 2"}{print $3}'`
endyyyy=`echo $stopdate | awk 'BEGIN{FIELDWIDTHS="4 2 2"}{print $1}'`
endmm=`echo $stopdate | awk 'BEGIN{FIELDWIDTHS="4 2 2"}{print $2}'`
enddd=`echo $stopdate | awk 'BEGIN{FIELDWIDTHS="4 2 2"}{print $3}'`

polygon=`awk 'BEGIN{ORS=","}{if($1!=">"){print $1, $2}}' $polyfile | sed 's/,$//' `

apiroot="https://scihub.copernicus.eu/dhus/search?"

# Always search IW mode:
querytags1=" AND sensoroperationalmode:IW"

# Check for placeholder arguements "-"
if [ $orbitdir=='-' ]; then
  querytags1="${querytags1}"
else
  querytags1="${querytags1} AND orbitdirection:$orbitdir"
fi

if [ $track == '-' ]; then
  querytags1="${querytags1}"
else
  querytags1="${querytags1} AND relativeorbitnumber:$track "
fi

querytags2="AND ingestiondate:[${startyyyy}-${startmm}-${startdd}T00:00:00.000Z TO ${endyyyy}-${endmm}-${enddd}T00:00:00.000Z] "

rm tmp1.txt tmp2.txt tmp3.txt
touch tmp1.txt tmp2.txt tmp3.txt
i=1
flag=1
while [ $flag == 1 ]; do
  rstart=`echo $i | awk '{print ($1-1)*100}'`
  querytags0="AND producttype:SLC&rows=${nmax}&start=${rstart}"
  wget --output-document=query_results.xml "${apiroot}q=footprint:\"Intersects(POLYGON((${polygon})))\" ${querytags1} ${querytags2} ${querytags0}"
  
  # Parse the xml output:
  grep -e 'filename' query_results.xml | sed 's/[\<,\>]/ /g' | awk '{print $3}' | sed 's/\.SAFE/\.zip/g' >> tmp1.txt
  grep "<link href=" query_results.xml | sed 's/\"/ /g' | awk '{print $3}' >> tmp2.txt
  grep 'name="relativeorbitnumber' query_results.xml | sed 's/[\<,\>]/ /g' | awk '{print $3}' >> tmp3.txt

  let i=i+1
  echo $i
  # Break the loop if the last query returned no useful rows
  if [ -z `grep -e 'filename' query_results.xml | awk '(NR==1){print 1}'` ]; then
    break
  fi
done
# wget "https://scihub.copernicus.eu/dhus/search?q=footprint:\"Intersects(POLYGON((-4.53 29.85,26.75 29.85,26.75 46.80,-4.53 46.80,-4.53 29.85)))\" "

paste tmp1.txt tmp2.txt tmp3.txt > $ziplist

rm tmp1.txt tmp2.txt tmp3.txt

awk '{x++}END{print x " zip files found."}' $ziplist




