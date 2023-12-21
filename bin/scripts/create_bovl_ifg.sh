#!/bin/bash

# 2023 : ML (with great support of Muhammet Nergizci)

if [ -z $1 ]; then 
 echo "USAGE: provide pair, and keep being in the frame folder"
 echo "e.g. create_bovl_ifg.sh 20230115_20230127 [1] ... to generate bovl ddiff ifg between those dates. Optional param [1] would filter the bovl"
 exit
fi

filterit=0
if [ ! -z $2 ]; then
 if [ $2 == 1 ]; then filterit=1;
 else
  echo "provided unexpected extra parameter, cancelling now"; exit
 fi
fi

source $LiCSARpath/lib/LiCSAR_bash_lib.sh

m=`echo $1 | cut -d '_' -f1`
s=`echo $1 | cut -d '_' -f2`
pair=$1
master=`get_master`
frame=`pwd`
frame=`basename $frame`
if [ `echo $frame | wc -m` != 18 ]; then echo "not in frame folder"; exit; fi

## first of all, do tabs:
if [ ! -f tab/$m'R_tab' ]; then
createSLCtab_frame RSLC/$m/$m rslc $frame > tab/$m'R_tab'
fi
if [ ! -f tab/$s'R_tab' ]; then
createSLCtab_frame RSLC/$s/$s rslc $frame > tab/$s'R_tab'
fi

mkdir -p IFG/$pair
stab=tab/$s'R_tab'
mtab=tab/$m'R_tab'
mastertab=tab/$master'R_tab'
if [ ! -f $mastertab ]; then
 createSLCtab_frame RSLC/$master/$master rslc $frame > $mastertab
fi

movlfile=IFG/$pair/$m'_overlaps'
sovlfile=IFG/$pair/$s'_overlaps'
# 1. extract burst overlaps, output *fwd.slc, fwd.slc.par, bwd.slc, bwd.slc.par
echo "extracting burst overlaps"
ScanSAR_burst_overlap $mtab $movlfile 20 4 0 0 $mastertab >/dev/null
ScanSAR_burst_overlap $stab $sovlfile 20 4 0 0 $mastertab >/dev/null


#3. Create offset, output offset files of bwd and fwd

create_offset $movlfile.fwd.slc.par $sovlfile.fwd.slc.par IFG/$pair/ovfwd.offset 1 20 4 0 >/dev/null
create_offset $movlfile.bwd.slc.par $sovlfile.bwd.slc.par IFG/$pair/ovbwd.offset 1 20 4 0 >/dev/null

#4. Make ifgs of backwards and forwards
echo "creating bwd and fwd ifgs"
SLC_intf $movlfile.bwd.slc $sovlfile.bwd.slc $movlfile.bwd.slc.par $sovlfile.bwd.slc.par IFG/$pair/ovbwd.offset IFG/$pair/bwd.ifg 20 4 0 - 0 0 >/dev/null
SLC_intf $movlfile.fwd.slc $sovlfile.fwd.slc $movlfile.fwd.slc.par $sovlfile.fwd.slc.par IFG/$pair/ovfwd.offset IFG/$pair/fwd.ifg 20 4 0 - 0 0 >/dev/null

#4.1 Make adf filtering to fwd adn bwd ifg
width=`get_value IFG/$pair/ovbwd.offset interferogram_width`
widthgeo=`get_value geo/EQA.dem_par width`

#if [ $filterit == 1 ]; then
# adf IFG/$pair/bwd.ifg IFG/$pair/bwd.adf.ifg IFG/$pair/bwd.adf.cc $width 1 - - - - - -
# adf IFG/$pair/fwd.ifg IFG/$pair/fwd.adf.ifg IFG/$pair/fwd.adf.cc $width 1 - - - - - -
 # to clean - not needed
# rm IFG/$pair/?wd.adf.cc
#fi
#5. Make their double difference

python3 -c "import numpy as np; a=np.fromfile('"IFG/$pair/"fwd.ifg', np.complex64).byteswap(); b=np.fromfile('"IFG/$pair/"bwd.ifg', np.complex64).byteswap(); c = a*np.conj(b);  c.byteswap().tofile('"IFG/$pair/"ddiff')"
# now, ddiff is after adf (!)
#python3 -c "import numpy as np; a=np.fromfile('"IFG/$pair/"fwd.adf.ifg', np.complex64).byteswap(); b=np.fromfile('"IFG/$pair/"bwd.adf.ifg', np.complex64).byteswap(); c = a*np.conj(b);  c.byteswap().tofile('"IFG/$pair/"ddiff')"


if [ $filterit == 1 ]; then
 echo "filtering"
 #5.2 Make adf double difference 
 #adf IFG/$pair/ddiff.adf IFG/$pair/ddiff.adf.adf IFG/$pair/ddiff.adf.adf.cc $width 1 - - - - - -
 adf IFG/$pair/ddiff IFG/$pair/ddiff.adf IFG/$pair/ddiff_coh_adf $width 1 - - - - - -  >/dev/null
fi
 #6. To create coherence
 # width=3447 
cc_wave IFG/$pair/ddiff RSLC/$m/$m.rslc.mli RSLC/$s/$s.rslc.mli IFG/$pair/ddiff_coh $width  >/dev/null


#7. create geotiffs
#width=`get_value RSLC/$master/$master.rslc.mli.par range_samples`
#width=`get_value IFG/$pair/ovbwd.offset interferogram_width`
#widthgeo=`get_value geo/EQA.dem_par width`

#
if [ $filterit == 1 ]; then
 cpx_to_real IFG/$pair/ddiff.adf IFG/$pair/ddiff_pha_adf $width 4  >/dev/null
fi
cpx_to_real IFG/$pair/ddiff IFG/$pair/ddiff_pha $width 4  >/dev/null

mkdir -p GEOC/$pair
echo "geocoding"
if [ $filterit == 1 ]; then
 geocode_back IFG/$pair/ddiff_coh_adf $width geo/$master.lt_fine GEOC/$pair/ddiff_coh_adf.geo $widthgeo - 0 0  >/dev/null
 geocode_back IFG/$pair/ddiff_pha_adf $width geo/$master.lt_fine GEOC/$pair/ddiff_pha_adf.geo $widthgeo - 0 0  >/dev/null
fi
geocode_back IFG/$pair/ddiff_coh $width geo/$master.lt_fine GEOC/$pair/ddiff_coh.geo $widthgeo - 0 0  >/dev/null
geocode_back IFG/$pair/ddiff_pha $width geo/$master.lt_fine GEOC/$pair/ddiff_pha.geo $widthgeo - 0 0  >/dev/null
#geocode_back IFG/$pair/ddiff_pha.adf $width geo/$master.lt_fine IFG/$pair/ddiff_pha.adf.geo $widthgeo - 0 0
if [ $filterit == 1 ]; then
 data2geotiff geo/EQA.dem_par GEOC/$pair/ddiff_coh_adf.geo 2 GEOC/$pair/$pair.geo.bovldiff.adf.cc.tif 0.0  >/dev/null
 data2geotiff geo/EQA.dem_par GEOC/$pair/ddiff_pha_adf.geo 2 GEOC/$pair/$pair.geo.bovldiff.adf.tif 0.0  >/dev/null
fi
data2geotiff geo/EQA.dem_par GEOC/$pair/ddiff_coh.geo 2 GEOC/$pair/$pair.geo.bovldiff.cc.tif 0.0  >/dev/null
data2geotiff geo/EQA.dem_par GEOC/$pair/ddiff_pha.geo 2 GEOC/$pair/$pair.geo.bovldiff.tif 0.0  >/dev/null
#data2geotiff geo/EQA.dem_par IFG/$pair/ddiff_pha.adf.geo 2 GEOC/$pair/$pair.geo.bovldiff.tif 0.0

## cc
#geocode_back IFG/$pair/ddiff.adf.cc $width geo/$master.lt_fine IFG/$pair/ddiff.adf.cc.geo $widthgeo - 0 0
#data2geotiff geo/EQA.dem_par IFG/$pair/ddiff.adf.cc.geo 2 GEOC/$pair/$pair.bovldiff.adf.cc.geo.tif 0.0


#8. create preview
if [ $filterit == 1 ]; then
 create_preview_wrapped GEOC/$pair/$pair.geo.bovldiff.adf.tif
else
 create_preview_wrapped GEOC/$pair/$pair.geo.bovldiff.tif
fi

#9. transform to mm
bovls_rad2mm.py $frame $pair

# 10. compress the final outputs
if [ $filterit == 1 ]; then tadf='adf.'; else tadf=''; fi
if [ -f GEOC/$pair/$pair.geo.bovldiff.$tadf'tif' ]; then
 mv GEOC/$pair/$pair.geo.bovldiff.$tadf'tif' GEOC/$pair/$pair.geo.bovldiff.$tadf'temp.tif';
 gdal_translate -co COMPRESS=DEFLATE -co PREDICTOR=3 -of GTiff -a_srs EPSG:4326 GEOC/$pair/$pair.geo.bovldiff.$tadf'temp.tif' GEOC/$pair/$pair.geo.bovldiff.$tadf'tif'
 rm GEOC/$pair/$pair.geo.bovldiff.$tadf'temp.tif'
fi
if [ -f GEOC/$pair/$pair.geo.bovldiff.$tadf'cc.tif' ]; then
 mv GEOC/$pair/$pair.geo.bovldiff.$tadf'cc.tif' GEOC/$pair/$pair.geo.bovldiff.$tadf'cc.temp.tif';
 gdal_translate -ot Byte -scale 0 1 0 255 -co COMPRESS=DEFLATE -co PREDICTOR=2 -of GTiff -a_srs EPSG:4326 GEOC/$pair/$pair.geo.bovldiff.$tadf'cc.temp.tif' GEOC/$pair/$pair.geo.bovldiff.$tadf'cc.tif'
 rm GEOC/$pair/$pair.geo.bovldiff.$tadf'cc.temp.tif'
fi
#20181029_20190103.geo.bovldiff.adf.mm.tif
if [ -f GEOC/$pair/$pair.geo.bovldiff.$tadf'mm.tif' ]; then
 mv GEOC/$pair/$pair.geo.bovldiff.$tadf'mm.tif' GEOC/$pair/$pair.geo.bovldiff.$tadf'mm.temp.tif';
 gdal_translate -co COMPRESS=DEFLATE -co PREDICTOR=3 -of GTiff -a_srs EPSG:4326 GEOC/$pair/$pair.geo.bovldiff.$tadf'mm.temp.tif' GEOC/$pair/$pair.geo.bovldiff.$tadf'mm.tif'
 rm GEOC/$pair/$pair.geo.bovldiff.$tadf'mm.temp.tif'
fi

# clean:
rm GEOC/$pair/ddiff* 2>/dev/null
echo "done - check preview of:"
ls GEOC/$pair/$pair.*bovldiff*.png
chmod 777 GEOC/$pair/*
