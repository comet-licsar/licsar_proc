#!/bin/bash
# 2019 - created by Daniel Juncu, edited by Milan Lazecky

if [ -z $3 ]; then 
 echo "This script will generate ifg from S1 stripmap images"
 echo "Parameters are: S1_STRIPMAP_ZIP_file_1 S1_STRIPMAP_ZIP_file_2 PARAMETER_file"
 echo "(absolute path for zip-files)"
 echo "(the processing will be executed in the current folder)"
 exit
fi

# source parameter file
source $3

echo "Extracting your stripmaps"
mkdir SLC; cd SLC
7za x $1
7za x $2
IF1DIR=`pwd`/`basename $1 .zip`.SAFE
IF2DIR=`pwd`/`basename $2 .zip`.SAFE
date1=`basename $1 .zip | cut -c 18-25`
date2=`basename $2 .zip | cut -c 18-25`
cd ..

mkdir $date1
cd $date1

par_S1_SLC ${IF1DIR}/measurement/*vv*.tiff ${IF1DIR}/annotation/*vv*.xml ${IF1DIR}/annotation/calibration/calibration*vv*.xml ${IF1DIR}/annotation/calibration/noise*vv*.xml ${date1}.slc.par ${date1}.slc
 
cd ..
mkdir $date2
cd $date2

par_S1_SLC ${IF2DIR}/measurement/*vv*.tiff ${IF2DIR}/annotation/*vv*.xml ${IF2DIR}/annotation/calibration/calibration*vv*.xml ${IF2DIR}/annotation/calibration/noise*vv*.xml ${date2}.slc.par ${date2}.slc

cd ..

#### pre-processing steps

### Manipulation of orbital state vectors
## generation of additional state vectors
# ORB_prop_SLC
## Modification of state vectors
# S1_OPOD_vec (<- problem: needs specific state vector file)

## Calibration: maybe not necessary for stripmap?


########
## generate offset parameter file
echo "======================================"
echo "ESTIMATE OFFSETS"
echo "======================================"
create_offset ${date1}/${date1}.slc.par ${date2}/${date2}.slc.par $date1'_'$date2.off 1 1 1 0


## estimation of offset using orbit information
init_offset_orbit ${date1}/${date1}.slc.par ${date2}/${date2}.slc.par $date1'_'$date2.off
## improve estimate (multi-looked)
init_offset ${date1}/${date1}.slc ${date2}/${date2}.slc ${date1}/${date1}.slc.par ${date2}/${date2}.slc.par $date1'_'$date2.off $rlks $azlks
## improve estimate (single-look)
init_offset ${date1}/${date1}.slc ${date2}/${date2}.slc ${date1}/${date1}.slc.par ${date2}/${date2}.slc.par $date1'_'$date2.off 1 1

## estimate field of offsets (vary parameters to improve results)
offset_pwr ${date1}/${date1}.slc ${date2}/${date2}.slc ${date1}/${date1}.slc.par ${date2}/${date2}.slc.par $date1'_'$date2.off offs ccp $offpwr_rwin $offpwr_azwin offsets $offpwr_n_ovr $offpwr_nr $offpwr_naz $offpwr_thres
## -OR-
# offset_SLC slc1 slc2 slc1.par slc2.par .off offs ccp 128 128 offsets 1 8 8 0.1

## compute offset polynomial
offset_fit offs ccp $date1'_'$date2.off coffs coffsets $offfit_coffs $offfit_coffsets $offfit_thres

## resample
echo "======================================"
echo "RESAMPLE"
echo "======================================"
SLC_interp ${date2}/${date2}.slc ${date1}/${date1}.slc.par ${date2}/${date2}.slc.par $date1'_'$date2.off ${date2}/${date2}.rslc ${date2}/${date2}.rslc.par 

## generate interferogram
echo "======================================"
echo "GENERATE INTERFEROGRAM"
echo "======================================"
SLC_intf ${date1}/${date1}.slc ${date2}/${date2}.rslc ${date1}/${date1}.slc.par ${date2}/${date2}.rslc.par $date1'_'$date2.off $date1'_'$date2.int $rlks $azlks $ifg_loff $ifg_nlines $ifg_sps_flg $ifg_azf_flg

## estimate baseline
echo "======================================"
echo "ESTIMATE BASELINE"
echo "======================================"
base_init ${date1}/${date1}.slc.par ${date2}/${date2}.slc.par $date1'_'$date2.off $date1'_'$date2.int $date1'_'$date2.base $base_mflag $base_nrfft $base_nazfft

## interferogram flattening
echo "======================================"
echo "FLATTEN INTERFEROGRAM"
echo "======================================"
ph_slope_base $date1'_'$date2.int ${date1}/${date1}.slc.par $date1'_'$date2.off $date1'_'$date2.base $date1'_'$date2.flt

## multi look images
echo "======================================"
echo "MULTI-LOOK SLC'S"
echo "======================================"
multi_look ${date1}/${date1}.slc ${date1}/${date1}.slc.par ${date1}/${date1}.mli ${date1}/${date1}.mli.par $rlks $azlks
multi_look ${date2}/${date2}.rslc ${date2}/${date2}.rslc.par ${date2}/${date2}.rmli ${date2}/${date2}.rmli.par $rlks $azlks

## estimate coherence
echo "======================================"
echo "ESTIMATE COHERENCE"
echo "======================================"
cc_wave $date1'_'$date2.flt ${date1}/${date1}.mli ${date2}/${date2}.rmli $date1'_'$date2.cc $cc_width $cc_bx $cc_by $cc_wflg

## adaptive filtering
# adapt_filt .int (or.flt?) .sm <width> # width: number of samples/row
## or 
echo "======================================"
echo "ADAPTIVE FILTERING"
echo "======================================"
adf $date1'_'$date2.flt $date1'_'$date2.flt_sm $date1'_'$date2.smcc $adf_width $adf_alpha $adf_nfft $adf_cc_win $adf_step $adf_loff $adf_nlines $adf_wfrac
