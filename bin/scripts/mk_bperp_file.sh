#!/bin/bash
# Make bperp_file for single master

## Master date
master=`ls geo/????????.hgt | cut -c 5-12`

## RSLC dates
ls -d RSLC/20??[01]?[0-3]?* | cut -c 6-13 | sort | uniq > rslc.list

## temporary working dir
tmpdir="RSLCtmp"
rm -rf $tmpdir
mkdir -p $tmpdir/$master

rm -f SLC_mosaics_tab
rm baselines  2>/dev/null

#copy all master par files first
for mp in `ls SLC/$master/*.slc.par`; do
 cp $mp $tmpdir/$master/`basename $mp | sed 's/slc/rslc/'`
done

echo $master $master 0.0000 0 > baselines

while read day; do
if [ $day != $master ]; then
 echo $day
 #put the rslc.par files to the temp folder
 ## If 7zipped
 if [ -f RSLC/${day}.7z ]; then
  if [ `ls -l  RSLC/${day}.7z | gawk {'print $5'}` -lt 1000 ]; then echo "file "$day'.7z is erroneous, removing'; rm RSLC/$day'.7z';
  else
   7za x -o$tmpdir RSLC/${day}.7z $day/$day.*rslc.par >/dev/null 2>/dev/null;
  fi
 elif [ -d RSLC/$day ]; then
  mkdir $tmpdir/$day
  cp RSLC/$day/$day.*rslc.par $tmpdir/$day/.
 fi
 #some par files are empty - remove them
 for par in `ls $tmpdir/$day/$day.*rslc.par`; do
  if [ `wc -l $par | gawk {'print $1'}` -eq 0 ]; then
   rm $par
  fi
 done
 par=`ls $tmpdir/$day/$day.*rslc.par 2>/dev/null | tail -n1 2>/dev/null`
 if [ -z $par ]; then echo "no valid par file exists for "$day;
 else
  ext=`echo $par | cut -d '/' -f3 | cut -c 9-`
  masterpar=$tmpdir/$master/$master$ext
  #perform the computation
  bperp=`base_orbit $masterpar $par - | grep perpendicular | gawk {'print $5'}`
  if [ ! -z $bperp ]; then
   btemp=`datediff $master $day`
   if [ $day -lt $master ]; then btemp="-"$btemp; fi
   echo $master" "$day" "$bperp" "$btemp >> baselines
  else
   echo "seems as erroneous rslc: "`pwd`/RSLC/$day
   echo "to debug:"
   echo "base_orbit $masterpar $par -"
   par2=`ls $tmpdir/$day/$day.*rslc.par 2>/dev/null | head -n1 2>/dev/null`
   if [ $par != $par2 ]; then
    #just last attempt to get the baselines -- already was working once!
    ext=`echo $par | cut -d '/' -f3 | cut -c 9-`
    masterpar=$tmpdir/$master/$master$ext
    bperp=`base_orbit $masterpar $par - | grep perpendicular | gawk {'print $5'}`
    if [ ! -z $bperp ]; then
     btemp=`datediff $master $day`
     if [ $day -lt $master ]; then btemp="-"$btemp; fi
     echo $master" "$day" "$bperp" "$btemp >> baselines
    fi
   fi
   #base_orbit $masterpar $par -
  fi
 fi
fi
done < rslc.list

## Make bperp_file wrt single master
#base_calc SLC_mosaics_tab $tmpdir/${master}/${master}.rslc.par bperp_file_tmp itab 0 >/dev/null 2>/dev/null

## Remove intermediate files
rm -f rslc.list
rm -rf $tmpdir
