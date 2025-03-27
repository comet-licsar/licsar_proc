#!/bin/bash

if [ -z $1 ]; then
 echo "Usage: ahb_postproc.sh FINALGEOCDIR frameID [ovrflag]"
 echo "e.g.:  ahb_postproc.sh GEOCml10GACOSmask 155D_02611_050400 1"
 echo "This script will perform extra processing and copying of outputs after finished LiCSBAS processing" # licsar2licsbas.sh"
 echo "PLEASE run it inside your frame directory (where you have GEOCmlX and TS_GEOCmlX directories)"
 echo "frameID MUST BE PROVIDED - please use the WHOLE FRAME ID as in LiCSAR. If you have the $frame'_2' etc and you are in such named folder, it will auto-identify it"
 echo "the FINALGEOCDIR must be the last one, for which you have TS_... folder existing"
 echo "please add '1' (as ovrflag) if you HAVE UPDATED MASK (and processed with step 16) to OVERWRITE EXISTING in AHB folder"
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
ovr=0
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
if [ $ovr -gt 0 ]; then
 echo "... eurasia-fixing and vstd fixing"
 LiCSBAS_vel_plate_motion.py -t TS_$geocd -f $frame -o $ahbdir/$outframe.vel_filt.mskd.eurasia.geo.tif --vstd_fix >/dev/null 2>/dev/null
 cp TS_$geocd/results/vstd_scaled.tif $ahbdir/$outframe.vstd_scaled.geo.tif
else
  if [ ! -f $ahbdir/$outframe.vstd_scaled.geo.tif ]; then
    echo "... eurasia-fixing and vstd fixing"
    LiCSBAS_vel_plate_motion.py -t TS_$geocd -f $frame -o $ahbdir/$outframe.vel_filt.mskd.eurasia.geo.tif --vstd_fix >/dev/null 2>/dev/null
    cp TS_$geocd/results/vstd_scaled.tif $ahbdir/$outframe.vstd_scaled.geo.tif
  fi
fi
# regenerate every time (assuming changes only in mask?)
for prd in vel.mskd vel.filt vel.filt.mskd; do
   if [ $ovr -gt 0 ]; then
 echo "... geotiffing "$prd
 LiCSBAS_flt2geotiff.py -i TS_$geocd/results/$prd -p $eqapar -o $ahbdir/$outframe.$prd.geo.tif >/dev/null 2>/dev/null
 cp TS_$geocd/results/$prd.png $ahbdir/$outframe.$prd.png
   else
     outprd=$ahbdir/$outframe.$prd.geo.tif
    if [ ! -f $outprd ]; then
     echo "... geotiffing "$prd
     LiCSBAS_flt2geotiff.py -i TS_$geocd/results/$prd -p $eqapar -o $outprd >/dev/null 2>/dev/null
     cp TS_$geocd/results/$prd.png $ahbdir/$outframe.$prd.png
    fi
   fi
done
# regenerate only if not existing:
for prd in coh_avg vstd loop_ph_avg_abs n_unw; do
outprd=$ahbdir/$outframe.$prd.geo.tif
if [ ! -f $outprd ]; then
 echo "... geotiffing "$prd
 LiCSBAS_flt2geotiff.py -i TS_$geocd/results/$prd -p $eqapar -o $outprd >/dev/null 2>/dev/null
 cp TS_$geocd/results/$prd.png $ahbdir/$outframe.$prd.png
fi
done
# coh_avg_XX - some are already generated... maybe only 24 days??
prd=`ls TS_$geocd/results/coh_avg_* | head -n 1 | rev | cut -d '/' -f 1 | rev | cut -d '.' -f 1`
outprd=$ahbdir/$outframe.$prd.geo.tif
if [ ! -f $outprd ]; then
 echo "... geotiffing "$prd
 LiCSBAS_flt2geotiff.py -i TS_$geocd/results/$prd -p $eqapar -o $outprd >/dev/null 2>/dev/null
 cp TS_$geocd/results/$prd.png $ahbdir/$outframe.$prd.png
fi
#
echo "... checking and copying cum.h5 file, EQA par etc as needed"
if [ $ovr -gt 0 ]; then
cp TS_$geocd/mask_ts.png $ahbdir/$outframe'_mask_ts.png'
cp TS_$geocd/network/network13.png $ahbdir/$outframe'_network.png'
cp TS_$geocd/cum.h5  $ahbdir/$outframe.cum.h5
fi
if [ ! -f $ahbdir/$outframe.EQA.dem_par ]; then
 cp $eqapar $ahbdir/$outframe.EQA.dem_par
fi
if [ ! -f $ahbdir/$outframe.cum.h5 ]; then
 cp TS_$geocd/cum.h5 $ahbdir/$outframe.cum.h5
fi
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

##### how to properly merge:
for fr in ????_?????_??????; do
  echo $fr
  cd $fr
  if [ ! -f $fr.vel.geo.tif ]; then echo $fr >> ../noveltif; else
  for x in iono tide gacos; do
    if [ -f $fr.$x.vel.tif ]; then
      if [ ! -f $fr.$x.vel.ok.tif ]; then
    gmt grdmath $fr.$x.vel.tif 0 DENAN $fr.vel.geo.tif 0 DENAN 0 NAN ISFINITE MUL 0 NAN = $fr.$x.vel.ok.tif=gd:GTiff
    fi; fi;
  done
  fi
  cd /gws/nopw/j04/nceo_geohazards_vol1/public/LiCSAR_products/AHB
done
# this would nanify values from ionovel where vel is 0 or nan


## how to do in jobs (run for iono where needed)
for x in ????_?????_??????; do
  if [ ! -f $x/$x.iono.vel.tif ]; then echo $x; cd $x;
  bsub2slurm.sh  -q short-serial -M 16384 -o ionolotus.o -e ionolotus.e /gws/nopw/j04/nceo_geohazards_vol1/public/LiCSAR_products/AHB/Milan/doinslurm.sh $x.cum.h5 $x;
  cd /gws/nopw/j04/nceo_geohazards_vol1/public/LiCSAR_products/AHB;
  fi;
done

# fixing cumfilt from yma:
for fr in 153D_04699_131413 167D_05077_131313 175A_04990_131313; do cd $fr;
mv $fr.cum.h5 backup.h5; mv cum_filt.h5 $fr.cum.h5;
ionoveltif=$fr.iono.vel.tif
if [ ! -f $ionoveltif ]; then
 eqa=`ls *EQA.dem_par | head -n 1`
 eqaok=$fr.EQA.dem_par
 if [ $eqa != $eqaok ]; then mv $eqa $eqaok; fi
if [ -f $eqaok ]; then
 epochsdir=$LiCSAR_public/`track_from_frame $fr`/$fr/epochs
 cumh=`ls *cum.h5 | head -n 1`
 cumhok=$fr.cum.h5
 if [ $cumh != $cumhok ]; then mv $cumh $cumhok; fi
 if [ ! -f $cumhok ]; then echo $fr >> ../missingionoset.rights; else
  # permissions
  chmod 775 $cumhok 2>/dev/null
  if [ `ls -alh $cumhok | cut -c 6` != 'w' ]; then
   echo "changing ownership"
   gwschown -n --no-warn earmla $cumhok >/dev/null 2>/dev/null
   gwschown --no-warn earmla $cumhok
   chmod 775 $cumhok
  fi
  python3 -c "from lics_tstools import *; correct_cum_from_tifs('"$cumhok"', '"$epochsdir"', 'tide.geo.tif', 1000, directcorrect = False)"
  python3 -c "from lics_tstools import *; correct_cum_from_tifs('"$cumhok"', '"$epochsdir"', 'geo.iono.code.tif', 55.465/(4*np.pi), directcorrect = False)"
  if [ ! -f $fr.gacos.vel.tif ]; then
    python3 -c "from lics_tstools import *; correct_cum_from_tifs('"$cumhok"', '"$epochsdir"', 'sltd.geo.tif', -55.465/(4*np.pi), directcorrect = False, newcumname = 'gacos')"
    LiCSBAS_cum2vel.py --datavar gacos -i $cumhok -o gacos;
    LiCSBAS_flt2geotiff.py -i gacos.vel -p $eqaok -o $fr.gacos.vel.tif
  fi
  LiCSBAS_cum2vel.py --datavar iono -i $cumhok -o iono
  LiCSBAS_flt2geotiff.py -i iono.vel -p $eqaok -o $ionoveltif
  LiCSBAS_cum2vel.py --datavar tide -i $cumhok -o tide
  LiCSBAS_flt2geotiff.py -i tide.vel -p $eqaok -o $fr.tide.vel.tif

 fi
else
 echo $fr >> ../missingionoset.noeqa
fi
fi

mv $fr.cum.h5 cum_filt.h5; mv backup.h5 $fr.cum.h5;
cd ..
done

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



#!/usr/bin/env python

import os
from pwd import getpwuid
import lics_unwrap as lu
import pandas as pd

def find_owner(filename):
    return getpwuid(os.stat(filename).st_uid).pw_name


def fix_enus(fr, ml=10, onlyrights=False):
    tr=str(int(fr[:3]))
    outdir='/gws/nopw/j04/nceo_geohazards_vol1/public/LiCSAR_products/AHB'
    outveltif=os.path.join(outdir, fr, fr+'.vel_filt.mskd.eurasia.geo.tif')
    if not os.path.exists(outveltif):
        print(fr+' not ready')
        rc = os.system('echo '+fr+' >> enusredo3.txt')
    else:
        for tc in ['E','N','U','hgt']:
            infile = os.path.join(os.environ['LiCSAR_public'], tr, fr[:17], 'metadata', fr[:17]+'.geo.'+tc+'.tif')
            outfile=os.path.join(outdir, fr, fr+'.'+tc+'.geo.tif')
            if onlyrights:
                if find_owner(outfile) != 'earmla':
                    #rc = os.system('gwschown --no-warn earmla '+outfile)
                    cmd='gwschown --no-warn earmla '+outfile
                    print(cmd)
            else:
                rc = os.system('chmod 777 '+outfile)
            if not onlyrights:
                xrda = lu.get_ml_hgt(infile, ml)
                lu.export_xr2tif(xrda, outfile, lonlat = True, dogdal = True, refto = outveltif)
                rc = os.system('chmod 444 '+outfile)


#testf='/gws/nopw/j04/nceo_geohazards_vol1/public/LiCSAR_products/AHB/173A_04952_131313/173A_04952_131313.E.geo.tif'
#os.system('gwschown -n --no-warn earmla '+testf)
f='enusredo.txt'
a=pd.read_csv(f,header=None)
frames = list(a[0].values)
for fr in frames:
    #print(fr)
    fix_enus(fr, onlyrights=True)

    173A_04952_131313