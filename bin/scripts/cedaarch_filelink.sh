#!/bin/bash 
#################################################################

usage() { 
  echo " " 1>&2; 
  echo "Usage: $0 " 1>&2; 
  echo "   This is to create file links from /neodc if the particular file exists (and is the same)" 1>&2;
  echo "   The param is a PNG or TIF file that should be in the \$LiCSAR_public folder.." 1>&2;
  echo "   e.g. cedaarch_filelink.sh \$file" 1>&2;
  echo " " 1>&2;
  exit 1; 
}

if [ -z $1 ]; then
  usage
else
  filepubdir=`realpath $1`
fi


function comparefiles() {
  cmp $1 $2 >/dev/null && echo "identical" || echo "different"
}

# basic checks
if [ -L $filepubdir ]; then
  #if [ `readlink $filepubdir | grep -c neodc` -eq 1 ]; then
  echo "the file is a link - skipping" # (may improve for neodc check)
  exit 1
fi
if [ ! -f $filepubdir ]; then
  echo "the file does not exist, exiting"
  exit 1
fi
ext=`echo "${filepubdir##*.}"`  # from https://stackoverflow.com/questions/965053/extract-filename-and-extension-in-bash
if [[ ! "png tif" =~ $ext ]]; then
  echo "cheking only if the extension is tif or png - the file is "$ext;
  exit 1
fi

# construct expected /neodc path
# /gws/nopw/j04/nceo_geohazards_vol1/public/LiCSAR_products.public/41/041A_05430_131313/interferograms/20170724_20210820/20170724_20210820.geo.unw.tif
filrev=`echo $filepubdir | rev`
track=`echo $filrev | cut -d '/' -f 5 | rev`
frame=`echo $filrev | cut -d '/' -f 4 | rev`
ddir=`echo $filrev | cut -d '/' -f 3 | rev`
pairep=`echo $filrev | cut -d '/' -f 2 | rev`
filename=`basename $filepubdir`

if [ ! $ddir == interferograms ] && [ ! $ddir == epochs ]; then
  echo "the script will work only for epochs or interferograms directories - cancelling for now"
  exit 1
fi

if [ $track -lt 10 ]; then cedatr=0$track; else cedatr=$track; fi
neodcframepath=/neodc/comet/data/licsar_products/$cedatr/$frame
# on CEDA Archive, we have ifgs directly in the frame dir, but epochs (or metadata) are under this..
if [ $ddir != interferograms ]; then
  neodcframepath=$neodcframepath/$ddir
fi

fileneodc=$neodcframepath/$pairep/$filename

# finally create link if this is identical file to the one in /neodc
if [ -f $fileneodc ]; then
if [ `comparefiles $filepubdir $fileneodc` == 'identical' ]; then  # echo "linking from /neodc";
   rm -f $filepubdir;
   if [ ! -f $filepubdir ]; then
     ln -s $fileneodc $filepubdir;
     chmod 755 $filepubdir;
   else
     echo "WARNING, could not delete file" $filepubdir;
     exit 1
   fi
fi
fi

