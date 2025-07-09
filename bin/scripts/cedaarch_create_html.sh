#!/bin/bash 
#################################################################
# This is to create html with links to either $LiCSAR_web or CEDA Archive
#################################################################

usage() { 
  echo " " 1>&2; 
  echo "Usage: $0 " 1>&2; 
  echo "   This is to create html with links to either \$LiCSAR_public or CEDA Archive " 1>&2;
  echo "   e.g. cedaarch_create_html.sh \$frame \$pair" 1>&2;
  echo "   or cedaarch_create_html.sh \$frame \$epoch epochs" 1>&2;
  echo "   (third param is core directory, by default: interferograms)" 1>&2;
  echo " " 1>&2; 
  exit 1; 
}

if [ -z $2 ]; then
  usage
else
  frame=$1
  pairep=$2
fi

ddir=interferograms
if [ ! -z $3 ]; then ddir=$2; fi
tr=`track_from_frame $frame`

indir=$LiCSAR_public/$tr/$frame/$ddir/$pairep
if [ ! -d $indir ]; then
  echo "The source dir does not exist: "$indir
  echo "cancelling cedaarch_create_html.sh script"
  exit
fi

outdir=$LiCSAR_web/$tr/$frame
outhtml=$outdir/$ddir
if [ ! -d $outdir ]; then
  mkdir -p $LiCSAR_web/$tr/$frame
  cd $LiCSAR_web/$tr/$frame
  if [ -d ../../../LiCSAR_products/$tr/$frame/metadata ]; then
   ln -s ../../../LiCSAR_products/$tr/$frame/metadata
  fi
  cd -
fi

pubwebdir=
cedaarchwebdir="https://data.ceda.ac.uk/neodc/comet/data/licsar_products/"$tr/$frame

if [ $ddir != interferograms ]; then
  cedaarchwebdir=$cedaarchwebdir/$ddir
fi

for ext in png tif; do
  for f in `ls $indir/*.$ext`; do
    # if the file is symbolic link leading to neodc:

    # else just make href to $LiCSAR_public