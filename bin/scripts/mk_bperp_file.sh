#!/bin/bash
#
# Make bperp_file for single master
#

## Master date
master=`ls geo/????????.hgt | cut -c 5-12`

## RSLC dates
ls -d RSLC/20??[01]?[0-3]?* | cut -c 6-13 | sort | uniq > rslc.list

## temporary working dir
tmpdir="RSLCtmp"
rm -rf $tmpdir
mkdir $tmpdir

## for each RSLC date
rm -f SLC_mosaics_tab
while read day; do
  rslc="$day/${day}.rslc"
  par="$day/${day}.rslc.par"

  ## If not 7zipped
  if [ -f RSLC/$par ]; then
    #echo "Link $par to $tmpdir"
    mkdir $tmpdir/$day
    cd $tmpdir/$day
    ln -s ../../RSLC/$par ./
    cd ../../
    echo "$tmpdir/$day/${day}.rslc $tmpdir/$day/${day}.rslc.par" >> SLC_mosaics_tab

  ## If 7zipped
  elif [ -f RSLC/${day}.7z ]; then
    #echo "Extract rslc.par for $day"
    ## Extract only rslc.par file
    7za x '-i!????????/????????.rslc.par' -o$tmpdir RSLC/${day}.7z

    if [ -f $tmpdir/$par ] ;then ## If succeed
      echo "$tmpdir/$day/${day}.rslc $tmpdir/$day/${day}.rslc.par" >> SLC_mosaics_tab
    fi

  fi

done < rslc.list

## Make bperp_file wrt single master
base_calc SLC_mosaics_tab $tmpdir/${master}/${master}.rslc.par bperp_file_tmp itab 0 >/dev/null 2>/dev/null

## Add 8th and 9th column
awk '{printf "%4d %8d %8d %7.2f %6.1f %6.1f %6.1f %7.2f %7.2f\n",$1,$2,$3,$4,$5,0,$5,0,$4}' bperp_file_tmp > bperp_file_old
awk '{printf "%8d %8d %7.2f %6.1f\n",$2,$3,$4,$5}' bperp_file_tmp > bperp_file

## Remove intermediate files
rm -f itab rslc.list SLC_mosaics_tab base_calc.log base.out bperp_file_tmp
rm -rf $tmpdir
