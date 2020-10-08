#!/bin/bash
#module load licsar_framebatch
#source /gws/smf/j04/nceo_geohazards/software/condalics/load_condalics.rc
month=`date '+%Y/%m'`
for path in `ls /neodc/sentinel1*/data/IW/L1_SLC/IPF_v*/$path -d`; do
 for i in `seq 1 7`; do
  arch2DB.py -d $path/0$i
 done
done
