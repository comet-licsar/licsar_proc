#!/bin/bash 
#################################################################

usage() {
  echo " " 1>&2;
  echo "Usage: $0 " 1>&2;
  echo "   This is to create file links from /neodc for all files in the (public) frame folder" 1>&2;
  echo "   The param is the frame ID" 1>&2;
  echo "   e.g. cedaarch_framelink.sh \$frame" 1>&2;
  echo " " 1>&2;
  exit 1;
}

if [ -z $1 ]; then
  usage
else
  frame=$1
fi

frpubdir=$LiCSAR_public/`track_from_frame $frame`/$frame
if [ ! -d $frpubdir ]; then
  echo "The frame directory does not exist:"
  echo $frpubdir
  exit 1
fi

 cd $frpubdir
 for fdir in `ls interferograms/????????_???????? epochs/?????? -d`; do
   for ext in tif png; do
     for f in `find $fdir -name '*.'$ext`; do
        cedaarch_filelink.sh $f
     done
   done
 done

 cd -

 cedaarch_create_html_frame.sh $frame