#!/bin/bash

if [ -z $1 ]; then
 echo "Usage e.g.: volq_get_volclip_path.sh [-n volcname] [-v volcID] "
 echo "try e.g. -n entale "
 exit
fi

while getopts ":n:v:" option; do
 case "${option}" in
  n ) vid=`python3 -c "import volcdb; volcid=int(volcdb.find_volcano_by_name('"$OPTARG"').volc_id); print(volcdb.get_volclip_vids(volcid)[0])" | tail -n 1`;
     ;;
  v ) vid=`python3 -c "import volcdb; print(volcdb.get_volclip_vids("$OPTARG")[0])" | tail -n 1`;
     ;;
 esac
done
shift $((OPTIND -1))

if [ -z $vid ]; then
  echo "volcano clip ID not found"
  exit
fi
echo $LiCSAR_volc/$vid