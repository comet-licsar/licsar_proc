

# script to reunwrap all and save as geo.unw_GACOS.tif
for aam in /gws/nopw/j04/nceo_geohazards_vol1/projects/LiCS/volc-proc/current/central_america /gws/nopw/j04/nceo_geohazards_vol1/projects/LiCS/volc-proc/current/south_america; do
cd $aam
pwd
for x in *_*_*_*; do
  echo $x
 frame=`echo $x | rev | cut -c -17 | rev`
 ddir=/work/scratch-pw3/licsar/earmla/batchdir/reunw/$x
 mkdir $ddir
 echo "from lics_unwrap import process_ifg_pair" > $ddir/todo.py
 echo "import os" >> $ddir/todo.py
 allprocessed=1
 for phatif in `ls $x/tif/*.diff_pha.tif `; do
   pp=`echo $phatif | sed 's/diff_pha/unw/'`
   procit=0
   if [ ! -f $pp.backup ]; then
     procit=1
     mv $pp $pp.backup
   elif [ ! -f $pp ]; then
     procit=1
   fi
   # one more check - if unw exists but is bad
   if [ $procit == 0 ]; then
     if [ `gdalinfo $pp | grep -c Description` == 0 ]; then rm $pp; procit=1; fi #
   fi
   if [ $procit == 1 ]; then
    allprocessed=0
    pair=`echo $pp | rev | cut -c 13-29 | rev`
    cohtif=`echo $phatif | sed 's/diff_pha/cc/'`
    echo "try:" >> $ddir/todo.py
    echo "    process_ifg_pair('"$phatif"', '"$cohtif"', procpairdir='/work/scratch-pw3/licsar/earmla/batchdir/reunw/"$x"', ml=1, fillby='nearest', thres=0.35, specmag=True, pre_detrend=False, outtif='"$pp"')" >> $ddir/todo.py
    echo "except:" >> $ddir/todo.py
    echo "    print('error in pair "$pair"')" >> $ddir/todo.py
    echo "" >> $ddir/todo.py
    echo "" >> $ddir/todo.py
   fi
 done
 chmod 777 $ddir/todo.py
 #python3 $ddir/todo.py # send to JASMIN!
 if [ $allprocessed == 0 ]; then
  bsub2slurm.sh -o $ddir/reunw.out -e $ddir/reunw.err -J reunw_$x -q short-serial -n 1 -W 12:59 -M 8200 python3 $ddir/todo.py # send to JASMIN!
 else
   echo "all unws already reprocessed, doing nothing"
 fi
done
done









# below by RR:



maule - update hires frames (gold)

VOLC_REGION=
VOLC_NAME=
VOLC_FRAME=`pwd`
VOLC_FRAME=`basename $VOLC_FRAME`
WEB_RESOLS='50 75 100'
OUT_DIR='/gws/nopw/j04/nceo_geohazards_vol1/projects/LiCS/volc-portal/processing_output/licsbas_gacos'

JSON_SCRIPT1='/gws/nopw/j04/nceo_geohazards_vol1/projects/LiCS/volc-portal/processing_scripts/licsbas/data_to_json.py'
# script for making JSON for web usage:
JSON_SCRIPT2='/gws/nopw/j04/nceo_geohazards_vol1/projects/LiCS/volc-portal/processing_scripts/licsbas/convert_json_for_web.py'


python ${JSON_SCRIPT1} */cum.h5 data.json
python ${JSON_SCRIPT1} */cum_filt.h5 data_filt.json


for WEB_RESOL in ${WEB_RESOLS}
do
  python ${JSON_SCRIPT2} data.json data_web_x${WEB_RESOL}.json ${WEB_RESOL}
  python ${JSON_SCRIPT2} data_filt.json data_web_x${WEB_RESOL}_filt.json ${WEB_RESOL}
done


\cp --preserve=timestamps data.json \
  ${OUT_DIR}/${VOLC_REGION}/${VOLC_NAME}_${VOLC_FRAME}.json
\cp --preserve=timestamps data_filt.json \
  ${OUT_DIR}/${VOLC_REGION}/${VOLC_NAME}_${VOLC_FRAME}_filt.json
for WEB_RESOL in ${WEB_RESOLS}
do
  \cp --preserve=timestamps data_web_x${WEB_RESOL}.json \
    ${OUT_DIR}/${VOLC_REGION}/${VOLC_NAME}_${VOLC_FRAME}_web_x${WEB_RESOL}.json
  \cp --preserve=timestamps  data_web_x${WEB_RESOL}_filt.json \
    ${OUT_DIR}/${VOLC_REGION}/${VOLC_NAME}_${VOLC_FRAME}_web_x${WEB_RESOL}_filt.json
done

\cp log/*.log \
  ${OUT_DIR}/${VOLC_REGION}/${VOLC_NAME}_${VOLC_FRAME}_latest.log


