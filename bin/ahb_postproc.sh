#!/bin/bash

if [ -z $1 ]; then
 echo "Usage: ahb_postproc.sh FINALGEOCDIR frameID"
 echo "e.g.:  ahb_postproc.sh GEOCml10GACOSmask 155D_02611_050400"
 echo "This script will perform extra processing and copying of outputs after finished LiCSBAS processing" # licsar2licsbas.sh"
 echo "PLEASE run it inside your frame directory (where you have GEOCmlX and TS_GEOCmlX directories)"
 echo "frameID MUST BE PROVIDED - please use the WHOLE FRAME ID as in LiCSAR. If you have the $frame'_2' etc and you are in such named folder, it will auto-identify it"
 echo "the FINALGEOCDIR must be the last one, for which you have TS_... folder existing"
 #echo "parameters: GEOCDIR frameID"
 #echo "e.g. GEOCml10GACOSmask 155D_02611_050400"
 echo "--------"
 echo "You should use this once you are happy with your final outputs (mask is good and you used step 16 after masking)."
 echo "NB: the command to run step 16 should be"
 echo ""
 echo "--------"
 echo "This script will:"
 echo " - regenerate vel wrt eurasia tif"
 echo " - create coh_avg, phase bias and other extra tifs"
 echo " - copy everything needed to the AHB directory"
 echo " - print the basic info"
 exit
fi

geocd=$1
frame=$2
outframe=`pwd | rev | cut -d '/' -f 1 | rev` # will solve the '_2' etc naming
if [ -z $frame ]; then echo "please provide the frame ID"; exit; fi
#frame=$outf`pwd | rev | cut -d '/' -f 1 | rev`; echo "assuming frame "$frame; fi
ahbdir=$LiCSAR_public/AHB/$outframe

if [ ! -d $ahbdir ]; then echo "ERROR, this directory does not exist: "$ahbdir; exit; fi

if [ ! -f TS_$geocd/cum_filt.h5 ]; then
  echo "Seems results do not exist here: TS_"$geocd
  echo "Cancelling"
  exit
fi

eqapar='TS_'$geocd/info/EQA.dem_par
if [ ! -f $eqapar ]; then echo "ERROR, this file does not exist - please fix:"; echo $eqapar; exit; fi

echo "Processing frame "$frame
echo "and storing to:"
echo $ahbdir
echo "----------"
echo "Regenerating all required geotiffs (and copying/renaming pngs) directly to the AHB folder"
#if [ ! -f $frame.vstd_scaled.geo.tif ]; then
 echo "... eurasia-fixing and vstd fixing"
 LiCSBAS_vel_plate_motion.py -t TS_$geocd -f $frame -o $ahbdir/$frame.vel_filt.mskd.eurasia.geo.tif --vstd_fix >/dev/null 2>/dev/null
 cp TS_$geocd/results/vstd_scaled.tif $ahbdir/$frame.vstd_scaled.geo.tif
# regenerate every time (assuming changes only in mask?)
for prd in vel.mskd vel_filt vel_filt.mskd; do
 echo "... geotiffing "$prd
 LiCSBAS_flt2geotiff.py -i TS_$geocd/results/$prd -p $eqapar -o $ahbdir/$frame.$prd.geo.tif >/dev/null 2>/dev/null
 cp TS_$geocd/results/$prd.png $ahbdir/$frame.$prd.png
done
# regenerate only if not existing:
for prd in coh_avg vstd loop_ph_avg_abs n_unw; do
outprd=$ahbdir/$frame.$prd.geo.tif
if [ ! -f $outprd ]; then
 echo "... geotiffing "$prd
 LiCSBAS_flt2geotiff.py -i TS_$geocd/results/$prd -p $eqapar -o $outprd >/dev/null 2>/dev/null
 cp TS_$geocd/results/$prd.png $ahbdir/$frame.$prd.png
fi
done
# coh_avg_XX - some are already generated... maybe only 24 days??
prd=`ls TS_$geocd/results/coh_avg_* | head -n 1 | rev | cut -d '/' -f 1 | rev | cut -d '.' -f 1`
outprd=$ahbdir/$frame.$prd.geo.tif
if [ ! -f $outprd ]; then
 echo "... geotiffing "$prd
 LiCSBAS_flt2geotiff.py -i TS_$geocd/results/$prd -p $eqapar -o $outprd >/dev/null 2>/dev/null
 cp TS_$geocd/results/$prd.png $ahbdir/$frame.$prd.png
fi
#
echo "... copying cum.h5 file, EQA par etc"
cp TS_$geocd/mask_ts.png $ahbdir/$frame'_mask_ts.png'
cp TS_$geocd/network/network13.png $ahbdir/$frame'_network.png'
cp TS_$geocd/cum.h5  $ahbdir/$frame.cum.h5
cp $eqapar $ahbdir/$frame.EQA.dem_par
echo ""
# nah, we will instead downsample the original ENUs and hgt + apply gdalwarp2match.py - I will do it (Milan)
#for x in U E N hgt; do
#  if [ ! -f $frame.$x.geo.tif ]; then
#    LiCSBAS_flt2geotiff.py -i $geocd/$x -p $geocd/EQA.dem_par -o $frame.$x.geo.tif
#  fi
#done

echo "done. Basic info on the dataset"
info=`ls TS_$geocd/info/13used_image.txt`
noim=`cat $info | wc -l`
frst=`head -n 1 $info `
lst=`tail -n 1 $info `
infoi=`ls TS_$geocd/info/13resid.txt`
noifgs=`cat $infoi | wc -l`
let noifgs=$noifgs-1
echo $frst'-'$lst','$noim','$noifgs
cat TS_$geocd/results/vstd_rescaling_parameters.txt

echo ""

exit

below are comments to the other things:

fr=''  # set the frame ID – this frame must be already inside the $LiCSAR_public/AHB directory
import lics_unwrap as lu
import glob, os
outdir='/gws/nopw/j04/nceo_geohazards_vol1/public/LiCSAR_products/AHB'
ml=10
tr=str(int(fr[:3]))
outveltif = os.path.join(outdir, fr, fr+'.vel_filt.mskd.eurasia.geo.tif')
checkerf = os.path.join(outdir, fr, fr+'.ENUS.ok')
if not os.path.exists(outveltif):
    print('not ready')
elif not os.path.exists(checkerf):
    os.system('touch '+checkerf)
    if os.path.exists(checkerf):
        for tc in ['E','N','U','hgt']:
            infile = os.path.join(os.environ['LiCSAR_public'], tr, fr, 'metadata', fr+'.geo.'+tc+'.tif')
            outfile=os.path.join(outdir, fr, fr+'.'+tc+'.geo.tif')
            xrda = lu.get_ml_hgt(infile, ml)
            lu.export_xr2tif(xrda, outfile, lonlat = True, dogdal = True, refto = outveltif)
