#!/bin/bash

usage() { 
  echo " " 1>&2; 
  echo "Usage: " 1>&2; 
  echo "  zips2cemszips.sh label_ziplist.full_list label_ziplist.list " 1>&2; 
  echo "   " 1>&2; 
  echo "  This program will convert a query zip list to a cems ziplist (with appropriate cems filepaths). " 1>&2; 
  echo "  " 1>&2; 
  exit 1; 
}

# Check if input file exists
if [ -f $1 ]; then 
  full_list=$1
else 
  echo "File not found: $1"
  usage 
fi

ziplist=$2

# Root directory at cems:
cemsroot1a=/neodc/sentinel1a/data/IW/L1_SLC/IPF_v2
cemsroot1b=/neodc/sentinel1b/data/IW/L1_SLC/IPF_v2

rm $ziplist; touch $ziplist;

echo "reading: $full_list"
while read listline; do
  echo $listline
  stub=`echo $listline | awk '{print $1}' | sed 's/\.zip//'`
  yyyy=${stub:17:4}
  mm=${stub:21:2}
  dd=${stub:23:2}
  sat=${stub:0:3}
  echo $stub $yyyy $mm $dd $sat
  if [ $sat == "S1A" ]; then cemsroot=$cemsroot1a; else cemsroot=$cemsroot1b; fi

  echo "${cemsroot}/${yyyy}/${mm}/${dd}/${stub}.zip" >> $ziplist
done < $full_list
