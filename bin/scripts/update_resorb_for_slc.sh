#!/bin/bash

if [ -z $1 ]; then echo "Parameter is filename, i.e. S1.....zip "; exit; fi

# e.g. S1A_IW_SLC__1SDV_20170301T214607_20170301T214634_015504_019780_9990.zip
FILENAME=$1
python3 -c "import orbit_lib as orb; orb.updateOrbForZipfile('"$FILENAME"')"

exit
# the below cannot work anymore after migrating to CDSE:

a=`echo $FILENAME | cut -c 18-25`
DATUM_START=`echo ${a:0:4}-${a:4:2}-${a:6:2}`
a=`echo $FILENAME | cut -c 27-32`
TIME_START=`echo ${a:0:2}:${a:2:2}:${a:4:2}`
a=`echo $FILENAME | cut -c 34-41`
DATUM_STOP=`echo ${a:0:4}-${a:4:2}-${a:6:2}`
a=`echo $FILENAME | cut -c 43-48`
TIME_STOP=`echo ${a:0:2}:${a:2:2}:${a:4:2}`

SAB=`echo $FILENAME | cut -d '_' -f1`

#     https://qc.sentinel1.eo.esa.int/api/v1/?product_type=AUX_RESORB
#echo https://qc.sentinel1.eo.esa.int/api/v1/?product_type=AUX_RESORB\&validity_stop__gt=2017-10-01T00:38:58&validity_start__lt=2017-10-01T00:39:25&ordering=-creation_date&page_size=10
#wget -O temp_resorb.txt https://qc.sentinel1.eo.esa.int/api/v1/?product_type=AUX_RESORB\&validity_stop__gt=${DATUM_STOP}T${TIME_STOP}\&validity_start__lt=${DATUM_START}T${TIME_START}\&ordering=-creation_date\&page_size=10 2>/dev/null
wget -O temp_resorb.txt https://qc.sentinel1.copernicus.eu/api/v1/?product_type=AUX_RESORB\&validity_stop__gt=${DATUM_STOP}T${TIME_STOP}\&validity_start__lt=${DATUM_START}T${TIME_START}\&ordering=-creation_date\&page_size=10 2>/dev/null

for x in `python -mjson.tool temp_resorb.txt | grep remote_url | cut -d '"' -f4`; do 
 wget -nc -O $ORB_DIR/$SAB/RESORB/`basename $x` $x 2>/dev/null
done

rm temp_resorb.txt
