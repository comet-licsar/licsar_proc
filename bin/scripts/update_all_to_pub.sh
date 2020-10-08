#!/bin/bash
#a fast correction - update all coherence maps and unw ifgs if they do not exist in public website..

#for additional coherence map workaround
module load doris

pubdir=/gws/nopw/j04/nceo_geohazards_vol1/public/LiCSAR_products
cat << EOF > ~/logs/tmp_procdirs.txt
/gws/nopw/j04/nceo_geohazards_vol2/LiCS/temp/volc
/gws/nopw/j04/nceo_geohazards_vol2/LiCS/temp/volc/frames/103
/gws/nopw/j04/nceo_geohazards_vol2/LiCS/temp/volc/frames/124
/gws/nopw/j04/nceo_geohazards_vol2/LiCS/temp/volc/frames/134
/gws/nopw/j04/nceo_geohazards_vol2/LiCS/temp/volc/frames/142
/gws/nopw/j04/nceo_geohazards_vol2/LiCS/temp/volc/frames/152
/gws/nopw/j04/nceo_geohazards_vol2/LiCS/temp/volc/frames/16
/gws/nopw/j04/nceo_geohazards_vol2/LiCS/temp/volc/frames/26
/gws/nopw/j04/nceo_geohazards_vol2/LiCS/temp/volc/frames/32
/gws/nopw/j04/nceo_geohazards_vol2/LiCS/temp/volc/frames/62
EOF

for fromdir in `cat ~/logs/tmp_procdirs.txt`; do

echo "Updating public website from "$fromdir
cd $fromdir
for frame in `ls [0-9]*_*_* -d`; do
frameno=`echo $frame | cut -c -3 | sed 's/^0//' | sed 's/^0//'`
echo "-------------"
echo "Frame "$frame
echo "-------------"
procdir=$fromdir/$frame
cd $procdir
if [ -f ${procdir}/geo/EQA.dem_par ]; then
if [ -d IFG ]; then
for ifg in `ls IFG/*/*.unw | cut -d '/' -f2`; do
#pokud existuje unw ifg v public, pak udelej jen coherence
if [ -f $pubdir/$frameno/$frame/interferograms/$ifg/$ifg.geo.unw.bmp ]; then
 echo "Updating only coherence map for "$ifg
 logfile=~/logs/coh_update.log
 master=`ls $procdir/geo/*[0-9].hgt | rev | cut -d '/' -f1 | rev | cut -d '.' -f1 | head -n1`
 width=`awk '$1 == "range_samples:" {print $2}' ${procdir}/RSLC/$master/$master.rslc.mli.par`;
 width_dem=`awk '$1 == "width:" {print $2}' ${procdir}/geo/EQA.dem_par`
 length_dem=`awk '$1 == "nlines:" {print $2}' ${procdir}/geo/EQA.dem_par`
 reducfac_dem=`echo $width_dem | awk '{if(int($1/2000) > 1) print int($1/2000); else print 1}'`
 geocode_back ${procdir}/IFG/${ifg}/${ifg}.cc $width ${procdir}/geo/$master.lt_fine ${procdir}/GEOC/${ifg}/${ifg}.geo.cc ${width_dem} ${length_dem} 1 0 >> $logfile
 data2geotiff ${procdir}/geo/EQA.dem_par ${procdir}/GEOC/${ifg}/${ifg}.geo.cc 2 ${procdir}/GEOC/${ifg}/${ifg}.geo.cc.tif 0.0  >> $logfile 2>/dev/null
 #rascc ${procdir}/GEOC/${ifg}/${ifg}.geo.cc - ${width_dem} - - - $reducfac_dem $reducfac_dem 0 1 - - - ${procdir}/GEOC/${ifg}/${ifg}.geo.cc_blk.bmp >> $logfile
 #cpxfiddle workaround to have grayscale coherence
 byteswap -o ${ifg}.geo.cc.tmp 4 ${procdir}/GEOC/${ifg}/${ifg}.geo.cc >/dev/null
 cpxfiddle -q normal -w ${width_dem} -f r4 -o sunraster -c gray -M $reducfac_dem/$reducfac_dem -r 0/0.9 ${ifg}.geo.cc.tmp > ${procdir}/GEOC/${ifg}/${ifg}.geo.cc_blk.ras 2>/dev/null
 convert ${procdir}/GEOC/${ifg}/${ifg}.geo.cc_blk.ras ${procdir}/GEOC/${ifg}/${ifg}.geo.cc_blk.bmp
 convert -resize 30% ${procdir}/GEOC/${ifg}/${ifg}.geo.cc_blk.bmp ${procdir}/GEOC/${ifg}/${ifg}.geo.cc.bmp
 rm -f ${procdir}/GEOC/${ifg}/${ifg}.geo.cc_blk.ras ${ifg}.geo.cc.tmp
 #moving to pubdir
  mkdir -p $pubdir/$frameno/$frame/interferograms/$ifg
 cp $procdir/GEOC/$ifg/$ifg.geo.cc.bmp $procdir/GEOC/$ifg/$ifg.geo.cc.tif $pubdir/$frameno/$frame/interferograms/$ifg/.
else
 #v opacnem pripade udelej vse
  echo "Generating whole files for "$ifg
  create_geoctiffs_to_pub.sh $procdir $ifg >/dev/null 2>/dev/null
  #moving to pubdir
  for toexp in cc.bmp cc.tif diff.bmp diff_mag.tif diff_pha.tif unw.bmp unw.tif disp.png; do
   mkdir -p $pubdir/$frameno/$frame/interferograms/$ifg
   cp $procdir/GEOC/$ifg/$ifg.geo.$toexp $pubdir/$frameno/$frame/interferograms/$ifg/.
  done
fi
done
else
 echo "There is no IFG folder here. Skipping"
fi
fi
done
done

rm ~/logs/tmp_procdirs.txt
