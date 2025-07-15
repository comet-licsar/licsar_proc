#!/bin/bash 
#################################################################
# create htmls for whole frame
#################################################################

usage() { 
  echo " " 1>&2; 
  echo "Usage: $0 " 1>&2; 
  echo "   This is to create html with links to either \$LiCSAR_public or CEDA Archive " 1>&2;
  echo "   e.g. cedaarch_create_html_frame.sh \$frame" 1>&2;
  echo " " 1>&2; 
  exit 1; 
}

if [ -z $1 ]; then
  usage
else
  frame=$1
fi

tr=`track_from_frame $frame`
framedir=$LiCSAR_public/$tr/$frame

outdir=$LiCSAR_web/$tr/$frame
if [ ! -d $outdir ]; then
  mkdir -p $LiCSAR_web/$tr/$frame
  chmod 755 $LiCSAR_web/$tr/$frame
  cd $LiCSAR_web/$tr/$frame
  # we want to keep local links to metadata
  if [ -d $LiCSAR_public/$tr/$frame/metadata ]; then
   ln -s $LiCSAR_public/$tr/$frame/metadata
   chmod 755 metadata
  else
    echo "WARNING, no metadata folder exists for the given frame"
  fi
  cd -
fi

for ddir in interferograms epochs; do
  echo "linking "$ddir
  for e in `ls $framedir/$ddir | grep ^20`; do
    cedaarch_create_html.sh $frame $e $ddir
  done
done

