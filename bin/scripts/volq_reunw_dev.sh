

# script to reunwrap all and save as geo.unw_GACOS.tif
cd /gws/nopw/j04/nceo_geohazards_vol1/projects/LiCS/volc-proc/current/central_america
for x in *_*_*_*; do
 frame=`echo $x | rev | cut -c -17 | rev`
 echo "from lics_unwrap import process_ifg_pair" > todo.py
 for pp in `ls $x/tif/*.unw.tif `; do
  pair=`echo $pp | rev | cut -c 13-29 | rev`
  mv $pp $pp.backup
  mkdir /work/scratch-pw3/licsar/earmla/batchdir/reunw/$x
  echo "process_ifg_pair('"$phatif"', '"$cohtif"', procpairdir='/work/scratch-pw3/licsar/earmla/batchdir/reunw/"$x"', ml=1, fillby='nearest', thres=0.35, specmag=True, pre_detrend=False, outtif='"$pp"')" >> todo.py



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


