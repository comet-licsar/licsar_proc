#!/bin/bash
outdir=/gws/pw/j05/comet_lics/backup_2022_IFG

cd $LiCSAR_procdir
for tr in `seq 1 175`; do
 for frame in `ls $tr`; do 
   if [ `ls $tr/$frame/IFG 2>/dev/null | wc -l` != 0 ]; then
    echo $frame
    mkdir -p $outdir/$tr/$frame/IFG
    mv $tr/$frame/IFG/* $outdir/$tr/$frame/IFG/.
   fi
 done
done
