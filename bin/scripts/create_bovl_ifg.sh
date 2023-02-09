#!/bin/bash

# 2023 : ML (with great support of Muhammet Nergizci)

if [ -z $1 ]; then 
 echo "USAGE: provide pair, and keep being in the frame folder"
 echo "e.g. licsar_offset_tracking_pair.sh 20230115_20230127 ... to do px offset tracking between those dates"
 exit
fi

source $LiCSARpath/lib/LiCSAR_bash_lib.sh

m=`echo $1 | cut -d '_' -f1`
s=`echo $1 | cut -d '_' -f2`
pair=$1
master=`get_master`
frame=`pwd`
frame=`basename $frame`

## first of all, do tabs:
if [ `ls RSLC/$master/$master.IW?.rslc | wc -l` -eq 3 ]; then
createSLCtab RSLC/$m/$m rslc 1 3 > tab/$m'R_tab'
createSLCtab RSLC/$s/$s rslc 1 3 > tab/$s'R_tab'
else
echo "for now, only full IWs are supported - check this script (and improve if needed)"
exit
fi


mkdir -p IFG/$pair
stab=tab/$s'R_tab'
mtab=tab/$m'R_tab'
mastertab=tab/$master'R_tab'
if [ ! -f $mastertab ]; then
 createSLCtab RSLC/$master/$master rslc 1 3 > $mastertab
fi

movlfile=IFG/$pair/$m'_overlaps'
sovlfile=IFG/$pair/$s'_overlaps'
# 1. extract burst overlaps, output *fwd.slc, fwd.slc.par, bwd.slc, bwd.slc.par
ScanSAR_burst_overlap $mtab $movlfile 20 4 0 0 $mastertab
ScanSAR_burst_overlap $stab $sovlfile 20 4 0 0 $mastertab


#3. Create offset, output offset files of bwd and fwd

create_offset $movlfile.fwd.slc.par $sovlfile.fwd.slc.par IFG/$pair/ovfwd.offset 1 20 4 0
create_offset $movlfile.bwd.slc.par $sovlfile.bwd.slc.par IFG/$pair/ovbwd.offset 1 20 4 0

#4. Make ifgs of backwards and forwards

SLC_intf $movlfile.bwd.slc $sovlfile.bwd.slc $movlfile.bwd.slc.par $sovlfile.bwd.slc.par IFG/$pair/ovbwd.offset IFG/$pair/bwd.ifg 20 4 0 - 0 0
SLC_intf $movlfile.fwd.slc $sovlfile.fwd.slc $movlfile.fwd.slc.par $sovlfile.fwd.slc.par IFG/$pair/ovfwd.offset IFG/$pair/fwd.ifg 20 4 0 - 0 0

#5. Make their double difference

python3 -c "import numpy as np; a=np.fromfile('"IFG/$pair/"fwd.ifg', np.complex64).byteswap(); b=np.fromfile('"IFG/$pair/"bwd.ifg', np.complex64).byteswap(); c = a*np.conj(b);  c.byteswap().tofile('"IFG/$pair/"ddiff')"


#6. To create coherence

# width=3447
# cc_wave ddiff RSLC/20230117/20230117.rslc.mli RSLC/20230129/20230129.rslc.mli ddiff_coh $width


#7. create geotiffs
#width=`get_value RSLC/$master/$master.rslc.mli.par range_samples`
width=`get_value IFG/$pair/ovbwd.offset interferogram_width`
widthgeo=`get_value geo/EQA.dem_par width`

cpx_to_real IFG/$pair/ddiff IFG/$pair/ddiff_pha $width 4

mkdir -p GEOC/$pair
#geocode_back ddiff_coh $width geo/20220919.lt_fine ddiff_coh.geo 3867 - 0 0
geocode_back IFG/$pair/ddiff_pha $width geo/$master.lt_fine GEOC/$pair/ddiff_pha.geo $widthgeo - 0 0
#data2geotiff geo/EQA.dem_par ddiff_coh.geo 2 ddiff_coh.geo.tif 0.0
data2geotiff geo/EQA.dem_par GEOC/$pair/ddiff_pha.geo 2 GEOC/$pair/$pair.bovldiff.geo.tif 0.0

#8. create preview
create_preview_wrapped GEOC/$pair/$pair.bovldiff.geo.tif

echo "done - check preview of:"
ls GEOC/$pair/$pair.bovldiff.geo.png

