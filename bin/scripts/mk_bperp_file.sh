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

parcode=`ls RSLC/$master/*rslc.par | head -n1 | cut -d '.' -f2`
#parcode is e.g. IW1

## for each RSLC date
rm -f SLC_mosaics_tab
while read day; do
  rslc="$day/${day}.rslc"
  par="$day/${day}.${parcode}.rslc.par"

  ## If 7zipped
  if [ -f RSLC/${day}.7z ]; then
    #echo "Extract rslc.par for $day"
    ## Extract only rslc.par file
    #par=`7za l RSLC/${day}.7z | grep rslc.par | head -n1 | rev | gawk {'print $1'} | rev`
    #7za x '-i!????????/????????.rslc.par' -o$tmpdir RSLC/${day}.7z >/dev/null 2>/dev/null
    7za x -o$tmpdir RSLC/${day}.7z $par >/dev/null 2>/dev/null
    #if [ -f $tmpdir/$par ] ;then ## If succeed
    #  echo "$tmpdir/$day/${day}.rslc $tmpdir/$day/${day}.rslc.par" >> SLC_mosaics_tab
    #fi

  ## If not 7zipped
  elif [ -f RSLC/$par ]; then
    #par=`ls $day/${day}*rslc.par | head -n1`
    #echo "Link $par to $tmpdir"
    mkdir $tmpdir/$day
    cd $tmpdir/$day
    ln -s ../../RSLC/$par ./
    cd ../../
    #echo "$tmpdir/$day/${day}.rslc $tmpdir/$day/${day}.rslc.par" >> SLC_mosaics_tab
  fi

done < rslc.list

#removing the old file..
#rm bperp_file 2>/dev/null
rm baselines  2>/dev/null

## Make baselines file wrt single master
for image in `cat rslc.list`; do
  bperp=`base_orbit $tmpdir/${master}/${master}.$parcode.rslc.par $tmpdir/$image/$image.$parcode.rslc.par - 2>/dev/null | grep perpendicular | gawk {'print $5'}`
 if [ ! -z $bperp ]; then
  btemp=`datediff $master $image`
  if [ $image -lt $master ]; then btemp="-"$btemp; fi
  echo $master" "$image" "$bperp" "$btemp >> baselines
 else
  echo "seems as erroneous rslc: "`pwd`/RSLC/$image
 fi
done

## Make bperp_file wrt single master
#base_calc SLC_mosaics_tab $tmpdir/${master}/${master}.rslc.par bperp_file_tmp itab 0 >/dev/null 2>/dev/null

 
## Add 8th and 9th column
#awk '{printf "%4d %8d %8d %7.2f %6.1f %6.1f %6.1f %7.2f %7.2f\n",$1,$2,$3,$4,$5,0,$5,0,$4}' bperp_file_tmp > bperp_file_old
#awk '{printf "%8d %8d %7.2f %6.1f\n",$2,$3,$4,$5}' bperp_file_tmp > bperp_file

## Remove intermediate files
#rm -f itab rslc.list SLC_mosaics_tab base_calc.log base.out bperp_file_tmp
rm rslc.list 
rm -rf $tmpdir
