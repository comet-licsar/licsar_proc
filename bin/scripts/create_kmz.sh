#!/bin/bash

if [ -z $1 ]; then echo "parameter is full path to the folder containing geotiffs";
echo "e.g. /gws/nopw/j04/nceo_geohazards_vol1/public/LiCSAR_products/58/058A_05279_131311/products/20190317_20190323"
exit;
fi

pair=`basename $1`
if [ ! `echo $pair | wc -m` -eq 18 ]; then 
 echo "wrong path - it should end on pairname, e.g. 20190101_20190107"
 exit
fi

if [ ! -f $1/$pair.geo.unw.bmp ] || [ ! -f $1/$pair.geo.diff.bmp ] || [ ! -f $1/$pair.geo.unw.tif ]; then
 echo "some of the geocoded product does not exist here"
 exit
fi

### let's make it
echo "generating kmz file for pair "$pair
cd $1
unw_tif=$pair.geo.unw.tif
unw_bmp=$pair.geo.unw.bmp
#ifg_tif=$pair.geo.diff_pha.tif
ifg_bmp=$pair.geo.diff.bmp

#generating previews
cp -r $LiCSARpath/misc/kmlfiles files
convert $unw_bmp files/`basename $unw_bmp .bmp`.png
convert $ifg_bmp files/`basename $ifg_bmp .bmp`.png

#getting templates and fit it
cp $LiCSARpath/misc/template.kml $pair.kml

#sed -i 's///' $pair.kml
#sed -i 's/TYEART/'`date +'%Y'`'/' $pair.kml
sed -i 's/TFIRSTDATET/'`echo $pair | cut -d '_' -f1`'/' $pair.kml
sed -i 's/TLASTDATET/'`echo $pair | cut -d '_' -f2`'/' $pair.kml
sed -i 's/TFILEUNWT/'`basename $unw_bmp .bmp`.png'/' $pair.kml
sed -i 's/TFILEIFGT/'`basename $ifg_bmp .bmp`.png'/' $pair.kml

#coordinates
gdalinfo -stats $unw_tif | grep Coordinates -A4 | tail -n+2 | cut -d '(' -f2 | cut -d ')' -f1 > tmp.coord
rm tmp.ewsn 2>/dev/null
for x in `cat tmp.coord | cut -d ',' -f1 | sort -n -u`; do echo $x >> tmp.ewsn; done
for x in `cat tmp.coord | cut -d ',' -f2 | sort -n -u`; do echo $x >> tmp.ewsn; done

sed -i 's/TEASTT/'`sed '1q;d' tmp.ewsn`'/' $pair.kml
sed -i 's/TWESTT/'`sed '2q;d' tmp.ewsn`'/' $pair.kml
sed -i 's/TSOUTHT/'`sed '3q;d' tmp.ewsn`'/' $pair.kml
sed -i 's/TNORTHT/'`sed '4q;d' tmp.ewsn`'/' $pair.kml

zip -r $pair.zip $pair.kml files >/dev/null
mv $pair.zip $pair.kmz

rm tmp.coord tmp.ewsn
rm -r files $pair.kml $unw_tif.aux.xml
if [ `echo $1 | grep -c LiCSAR_products` -gt 0 ]; then
 echo "done. Being in public folder, the file should be accessible at:"
 echo "http://gws-access.ceda.ac.uk/public/nceo_geohazards/"`echo $1 | cut -c 43-`"/"$pair".kmz"
fi

