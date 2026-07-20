#!/bin/bash

pedrofld=/gws/ssde/j25a/nceo_geohazards/vol1/projects/COMET/DEEPVolc_test_page/LiCSBAS
outdir=/gws/ssde/j25a/nceo_geohazards/vol2/LiCS/temp/insar_proc/earmla/volc_ncs

echo "assuming you work from "$pedrofld

for volcidf in `ls */volc_id.txt`; do
  volcfld=`dirname $volcidf`
  volcid=`cat $volcidf`
  for filtnc in `ls $volcfld/*/TS*/cum_filt.h5`; do
    frameid=`echo $filtnc | cut -d '/' -f 2`
    outfile=$outdir/$volcid'_'$frameid.filt.nc
    if [ ! -s $outfile ]; then
      LiCSBAS_out2nc.py -i $filtnc -o $outfile --xcube --compress
    fi
  done

  for unfiltnc in  `ls $volcfld/*/TS*/cum.h5`; do
    frameid=`echo $unfiltnc | cut -d '/' -f 2`
    outfile=$outdir/$volcid'_'$frameid.nc
    if [ ! -s $outfile ]; then
      LiCSBAS_out2nc.py -i $unfiltnc -o $outfile --xcube --compress
    fi
  done
done

exit


hdf=TS*/cum_filt.h5
LiCSBAS_out2nc.py -i $hdf -o temp.nc --cf
nccopy -4 -d 4 -s -u temp.nc out.64.nc
