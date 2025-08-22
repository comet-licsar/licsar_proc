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
if [ ! -z $3 ]; then ddir=$3; fi
tr=`track_from_frame $frame`

indir=$LiCSAR_public/$tr/$frame/$ddir/$pairep
if [ ! -d $indir ]; then
  echo "The source dir does not exist: "$indir
  echo "cancelling cedaarch_create_html.sh script"
  exit
fi

if [ ! $ddir == interferograms ] && [ ! $ddir == epochs ]; then
  echo "unsupported directory name provided - cancelling for now"
  exit
fi

outdir=$LiCSAR_web/$tr/$frame

if [ ! -d $outdir ]; then
  mkdir -p $LiCSAR_web/$tr/$frame
  chmod 775 $LiCSAR_web/$tr/$frame
  cd $LiCSAR_web/$tr/$frame
  # we want to keep local links to metadata
  if [ -d $LiCSAR_public/$tr/$frame/metadata ]; then
   ln -s $LiCSAR_public/$tr/$frame/metadata
   chmod 775 metadata
  else
    echo "WARNING, no metadata folder exists for the given frame"
  fi
  cd -
fi
if [ ! -d $outdir/$ddir ]; then
  mkdir $outdir/$ddir
  chmod 775 $outdir/$ddir
fi

outhtml=$outdir/$ddir/$pairep
pubwebdir="https://gws-access.jasmin.ac.uk/public/nceo_geohazards/LiCSAR_products.public/"$tr/$frame/$ddir

if [ $tr -lt 10 ]; then cedatr=0$tr; else cedatr=$tr; fi
cedaarchwebdir="https://data.ceda.ac.uk/neodc/comet/data/licsar_products/"$cedatr/$frame
# on CEDA Archive, we have ifgs directly in the frame dir, but epochs or metadata are under this..
if [ $ddir != interferograms ]; then
  cedaarchwebdir=$cedaarchwebdir/$ddir
fi

rm -f $outhtml
for ext in png tif kmz; do   # we store only png and tif files in /neodc but need also kmzs
  for f in `ls $indir/*.$ext 2>/dev/null`; do
  if [ -f $f ]; then    # just in case this is symlink that does not exist...
    # if the file is symbolic link leading to neodc:
    if [ -L $f ]; then
      if [ `readlink $f | grep -c neodc` -eq 1 ]; then
         hrefdir=$cedaarchwebdir
      else
         echo "WARNING - file "$f" is a link but not to /neodc. May not work"
         hrefdir=$pubwebdir
      fi
    # else just make href to $LiCSAR_public
    else
      hrefdir=$pubwebdir
    fi
    bf=`basename $f`
    echo "<a href='"$hrefdir"/"$pairep"/"$bf"'>"$bf"</a><br />" >> $outhtml
  fi
  done
done

chmod 775 $outhtml
