#!/bin/bash

usage() { 
  echo " " 1>&2; 
  echo "Usage: " 1>&2; 
  echo "  zips2cemszips.sh label_ziplist.full_list label_ziplist.list [mode]" 1>&2; 
  echo "   " 1>&2; 
  echo "  This program will convert a query zip list to a cems ziplist (with appropriate cems filepaths). " 1>&2; 
  echo "  " 1>&2; 
  echo "   optional mode parameter can be SM, it is IW by default" 1>&2; 
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

#optional third parameter - for non-IW data
if [ ! -z $3 ]; then
 mode=$3
else
 mode='IW'
fi
# Root directory at cems:
cemsroot1a=/neodc/sentinel1a/data/$mode/L1_SLC/IPF_v
cemsroot1b=/neodc/sentinel1b/data/$mode/L1_SLC/IPF_v
cemsroot1c=/neodc/sentinel1c/data/$mode/L1_SLC/IPF_v

rm $ziplist 2>/dev/null; touch $ziplist;

echo "reading: $full_list"
while read listline; do
  echo $listline
  stub=`echo $listline | awk '{print $1}' | sed 's/\.zip//'`
  yyyy=${stub:17:4}
  mm=${stub:21:2}
  dd=${stub:23:2}
  sat=${stub:0:3}
  echo $stub $yyyy $mm $dd $sat
  if [ $yyyy$mm$dd -gt 20190624 ]; then vers=3; else vers=2; fi
  if [ $sat == "S1A" ]; then cemsroot=$cemsroot1a$vers;
  elif [ $sat == "S1B" ]; then cemsroot=$cemsroot1b$vers;
      else cemsroot=$cemsroot1c$vers; fi

  echo "${cemsroot}/${yyyy}/${mm}/${dd}/${stub}.zip" >> $ziplist
done < $full_list
