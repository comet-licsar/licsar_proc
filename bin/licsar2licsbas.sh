#!/bin/bash
# parameters: frame [startdate] [enddate]
# e.g. 155D_02611_050400 20141001 20200205
  
if [ -z $1 ]; then
 echo "parameters: frame [startdate] [enddate]"
 echo "e.g. 155D_02611_050400 20141001 20200205"
 echo "parameters:"
 echo "-- Basic parameters --"
 echo "-M 10 .... this will do extra multilooking (in this example, 10x multilooking)"
 echo "-g ....... use GACOS if available for at least half of the epochs (and use only ifgs with both epochs having GACOS correction, other will be skipped)"
 echo "-S ....... (obsolete but still kept parameter) strict mode - in case of GACOS, use it only if available for ALL ifgs"
 echo "-G lon1/lon2/lat1/lat2  .... clip to this AOI"
 echo "-u ....... use the reunwrapping procedure (useful if multilooking..)"
 echo "-- Additional corrections --"
 echo "-i ....... perform ionosphere correction (using CODE GIM)"
 echo "-e ....... perform solid Earth tides correction"
 echo "-O ....... outputs of external corrections will be stored in the datacube itself (currently with -i, -e, testing -g)"
 echo "-- Control over reunwrapping (with -u) --"
 echo "-c ....... if the reunwrapping is to be performed, use cascade (might be better, especially when with shores)"
 echo "-l ....... if the reunwrapping is to be performed, would do Gaussian lowpass filter (should be safe unless in tricky areas as islands; good to use by default)"
 #echo "-m ....... OBSOLETE: keeping alway on (would be off only with -s) with reunwrapping with Goldstein filter on (by default), use coh based on spectral magnitude - recommended param, please use this by default"
 echo "-s ....... if the reunwrapping is to be performed, use Gaussian smooth filtering (this will turn off Goldstein filter)"
 echo "-m ....... use GAMMA ADF for filtering if Goldstein filter is selected (does not work together with -s)"
 echo "-t 0.2 .. change consistence threshold to 0.2 (should be default, quite good for ML10.. but now default is set to 0 !) during reunwrapping"
 echo "-H ....... this will use hgt to support the (re-)unwrapping"
 echo "-f ....... during reunwrapping, store also as GeoTIFFs"
 echo "-- Control over LiCSBAS processing --"
 echo "-T ....... use testing version of LiCSBAS"
 echo "-d ....... use the dev parameters for the testing version of LiCSBAS (currently: this will use --nopngs and --nullify, in future this may also add --singular)"
 echo "-W ....... use WLS for the inversion (coherence-based)"
 echo "-r ....... perform deramping (degree 1)"
 echo "-L ....... use linear hgt correlation slope correction"
 echo "-- Processing tweaks --"
 echo "-h 14 .... set your own number of processing hours (14 by default)"
 echo "-P ....... prioritise, i.e. use comet queue instead of short-serial"
 echo "-n 1 ..... number of processors (by default: 1, used also for reunwrapping)"
 echo "(other params, for admins etc.)"
 echo "(-R ....... prioritise through comet responder)"
 #echo "-----------------"
 echo "-I ....... use ICAMS - NOTE THIS WILL NOT OUTPUT AS ICAMS-CORRECTED DATASET - it will only generate and load ICAMS to the cum.h5 as icams layer (plus gacos layer to allow diff correction manually)"
 echo "-a ....... use amplitude stability to subset pixels. Testing. Might be useful in challenging areas such as jungles"
 echo "-w ....... will avoid using landmask"
 #echo "-A ....... use ampcoh"
 echo "-F ....... will force start from filtered ifgs (MAY NOT WORK IN LOCAL (HIRES) OR CASCADE? ANYWAY NOT RECOMMENDED - bit more loop errors as briefly tested)"
 echo "Note: you may want to check https://comet-licsar.github.io/licsar_proc/index.html#reunwrapping-existing-interferograms"
 #echo "note: in case you combine -G and -u, the result will be in clip folder without GACOS! (still not smoothly combined reunw->licsbas, todo!)"  # updated on 2022-04-07
 #echo "(note: if you do -M 1, it will go for reprocessing using the cascade/multiscale unwrap approach - in testing, please give feedback to Milan)"
 exit
fi

thres=0
echo "WARNING - 2024-08-30: indeed current noise estimation is removing edges such as phase jumps. Setting thres=0 by default to avoid this but then we do not really mask noise - expect very different results.."
echo "(of course, you can return this by manually setting -t 0.23 for example .. sorry if you already did and see this message)"
echo "And other update (2024-10-18): using avg abs phase bias == 1 rad for masking in step 15"
dolocal=0
dogacos=0
multi=1
run_jasmin=1
#dogacos=0
hgts=0
clip=0
reunw=0
use_coh_stab=0
strict=0
keep_coh_debug=1
#LB_version=LiCSBAS
cascade=0
smooth=0
gammadf=0
lowpass=0
wls=0
cometdev=0
#specmag=0
# 2023 - keep it ON
specmag=1
nproc=1
ampstab=0
ampcoh=0
filtifg=0
que='short-serial'
#LB_version=LiCSBAS  # orig Yu Morishita's version
LB_version=licsbas_comet  # COMET LiCSBAS (main branch)
#LB_version=LiCSBAS_testing
bovls=0
setides=0
iono=0
deramp=0
storeext2cube=0
prelb_backup=0
lotushours=0
icams=0
landmask=1
maskbias=1
hgtcorrlicsbas=0
outifs=0

discmd="$0 $@"
while getopts ":M:h:HucTsdbSClWgmaAiIeFfOPRrLwkG:t:n:" option; do
 case "${option}" in
  h) lotushours=${OPTARG};
     ;;
  M) multi=${OPTARG};
     #shift
     ;;
  b) bovls=1;
     echo "setting to process bovl data. Not ready yet (although LiCSBAS is.. with a workaround)"
     ;;
  i) iono=1;
     prelb_backup=1;
    ;;
  r) deramp=1;
    ;;
  L) hgtcorrlicsbas=1;
    ;;
  e) setides=1;
    prelb_backup=1;
    ;;
  f) outifs=1;
    ;;
  n) nproc=${OPTARG};
     #que='par-single'; # unless changed to comet queue
     ;;
  H) hgts=1;
     #shift
     ;;
  w) landmask=0;
     ;;
  I) icams=1;
    echo "Warning, ICAMS will be only loaded into the cum.h5 datacube. Please chat with earmla - you may want to do gacos-icams difference and correct cum TS yourself"
    storeext2cube=1;
    ;;
  #m) specmag=1;
  #   echo 'parameter -m is obsolete, the Goldstein is always on now';
  #   ;;
  O) storeext2cube=1;
     ;;
  u) reunw=1;
     #shift
     ;;
  c) cascade=1;
     ;;
  g) dogacos=1;
     echo "checking gacos";
     ;;
  d) cometdev=1;
     ;;
  C) use_coh_stab=1;
     #shift
     ;;
  s) smooth=1;
     #shift
     ;;
  m) gammadf=1;
    ;;
  F) filtifg=1;
    ;;
  P) que='comet';
     ;;
  R) que='comet_responder';
     ;;
  W) wls=1;
     ;;
  l) lowpass=1;
     #shift
     ;;
  S) strict=1;
     #shift
     ;;
  k) keep_coh_debug=0;
     #shift
     ;;
  a) ampstab=1;
     ;;
  A) ampcoh=1;
     ;;
  T) LB_version=licsbas_comet_dev; # COMET LiCSBAS (dev branch)
     #shift
     ;;
  G) aoi=${OPTARG};
     clip=1;
     #echo "warning - the clipping will affect only LiCSBAS for now, so in case of ML, the clip will be done only AFTER all reunwrapping"
     #shift
     ;;
  t) thres=${OPTARG};
     #clip=1;
     #echo "warning - the clipping will affect only LiCSBAS for now, so in case of ML, the clip will be done only AFTER all reunwrapping"
     #shift
     ;;
 esac
done
shift $((OPTIND -1))

if [ $nproc -gt 1 ]; then
 if [ $que == 'short-serial' ]; then
  que='par-single';
 fi 
fi

# what extension is used as input?
if [ $reunw -gt 0 ]; then
  if [ $filtifg -gt 0 ]; then
    extofproc='diff_pha'
  else
    extofproc='diff_unfiltered_pha'
  fi
elif [ $bovls -gt 0 ]; then
  if [ $prelb_backup -gt 0 ]; then
    echo "Error - you tried running longwave signal removal from burst overlap interferograms - not implemented in this tool"
    exit
  fi
else
  extofproc='unw'
fi

if [ -d GEOC ]; then
 echo "warning - GEOC folder detected. will use its contents for processing, rather than link from LiCSAR_public"
 dolocal=1;
fi

# frame info
if [ -f sourceframe ]; then frame=`cat sourceframe`; echo "setting frame from sourceframe file to "$frame;
    else frame=$1;
fi

if [ `echo $frame | grep -c '_'` -lt 1 ]; then
 echo "this is not a frame - check your input parameters please, yet continuing"
fi



# INIT
source $LiCSARpath/lib/LiCSAR_bash_lib.sh
thisdir=`pwd`
if [ $dolocal == 0 ]; then
 if [ ! `pwd | rev | cut -d '/' -f1 | rev` == $frame ]; then
  mkdir $frame
  cd $frame
 fi;
fi

defdates=1
startdate=20141001
enddate=`date +%Y%m%d`
if [ ! -z $2 ]; then
 defdates=0
 startdate=$2
fi
if [ ! -z $3 ]; then
 defdates=0
 enddate=$3
fi

track=`echo $frame | cut -c -3 | sed 's/^0//' | sed 's/^0//'`

indir=$LiCSAR_public/$track/$frame/interferograms
epochdir=$LiCSAR_public/$track/$frame/epochs

if [ ! -d $epochdir ]; then
  echo "warning, epochs directory does not exist for this frame"
  if [ $prelb_backup -gt 0 ]; then
    echo "and you requested longwave signal estimation. Please make sure you provided correct frame ID, cancelling for now"
    exit
  fi
fi
metadir=$LiCSAR_public/$track/$frame/metadata
workdir=`pwd`
source $metadir/metadata.txt #this will bring 'master=' info


 # setting for JASMIN processing
 if [ $lotushours -eq 0 ]; then
  hours=14
  if [ $multi -eq 1 ]; then
   hours=23
  fi
 else
  hours=$lotushours
  echo "setting processing time to "$hours":59"
  echo "(expecting you know what you are doing - i.e. non-comet queue has limit up to 23:59)"
 fi
 memm=8192 # requesting 8GB RAM per processor
 let memmfull=$nproc'*'$memm
 
 
mkdir GEOC 2>/dev/null
if [ $dogacos == 1 ]; then
  mkdir -p GACOS # 2>/dev/null
  #echo "debug gacos 1"
fi
cd GEOC
for meta in E N U hgt; do
 if [ -f lookangles/$master.geo.$meta.tif ]; then
  echo "getting metafiles from GEOC/lookangles - might need updating in future"
  ln -s `pwd`/lookangles/$master.geo.$meta.tif `pwd`/$master.geo.$meta.tif
 else
  ln -s $metadir/$frame.geo.$meta.tif
 fi
done
if [ $dolocal == 1 ]; then
 if [ ! -f $frame.geo.hgt.tif ]; then
  #echo "warning, doing local proc only - avoiding possible problems (mismatch frame vs local data) by just removing the linked ENU files"
  # 2024 - what problems? TODO (becoming forgetful..)
  rm $frame.geo.?.tif
 fi
fi
ln -s $metadir/baselines

if [ -f $workdir/GEOC.MLI/$master/$master.geo.mli.tif ]; then
 ln -s $workdir/GEOC.MLI/$master/$master.geo.mli.tif
else
if [ ! -f $epochdir/$master/$master.geo.mli.tif ]; then
 if [ -d $LiCSAR_procdir/$track/$frame ]; then
 echo "regenerating missing master MLI geotiff"
 create_geoctiffs_to_pub.sh -M $LiCSAR_procdir/$track/$frame $master
 mkdir -p $epochdir/$master
 mv $LiCSAR_procdir/$track/$frame/GEOC.MLI/$master/* $epochdir/$master/.
 rmdir $LiCSAR_procdir/$track/$frame/GEOC.MLI/$master
 else
  echo "warning, frame dir does not exist, not linking master mli file"
 fi
fi
if [ -f $epochdir/$master/$master.geo.mli.tif ]; then
 ln -s $epochdir/$master/$master.geo.mli.tif
else
 echo "warning, primary epoch not in public dir. trying to use other, expect possible size issues.."
 lastimg=`ls $epochdir/*/*.geo.mli.tif | tail -n1`
 ln -s $lastimg
fi
fi


if [ $dogacos == 1 ]; then
for epoch in `ls $epochdir`; do
  if [ $epoch -ge $startdate ] && [ $epoch -le $enddate ]; then
    gacosfile=$epochdir/$epoch/$epoch.sltd.geo.tif
    if [ -f $gacosfile ]; then
     ln -s $gacosfile $workdir/GACOS/$epoch.sltd.geo.tif 2>/dev/null
    fi
  fi
done
fi

if [ $dolocal == 0 ]; then
  disdir=`pwd`
  echo "Linking tif files from the LiCSAR_public directory"
  ls $indir | grep '_' > tmp.ifgs
if [ ! -z $2 ]; then
  echo "limiting the dataset to dates between "$startdate" and "$enddate
    #cp $epochdir/$epoch/$epoch.sltd.geo.tif ../GACOS/. 2>/dev/null
  for pair in `cat tmp.ifgs`; do
   if [ `echo $pair | cut -d '_' -f1` -ge $startdate ]; then
    if [ `echo $pair | cut -d '_' -f2` -le $enddate ]; then
      #ln -s $ifg;
      if [ ! -d $pair ]; then
        mkdir -p $pair
        cd $pair
        for ff in `ls $indir/$pair/*tif`; do
          ln -s $ff
        done
        cd $disdir
      fi
    fi
   fi
  done
else
 for pair in `cat tmp.ifgs`; do
   if [ ! -d $pair ]; then
    mkdir -p $pair
    cd $pair
    for ff in `ls $indir/$pair/*tif`; do
      ln -s $ff
    done
    cd $disdir
   fi
 done
fi
 # once done, switch to 'local' so we can use those linked data further
 rm tmp.ifgs
 dolocal=1
fi

# restore the backed up prev tifs if any
echo "checking on temporary backup (in case of another run in the same dir)"
noprelb=`ls 20*/*prelb.tif 2>/dev/null | wc -l`
if [ $noprelb -gt 0 ]; then
for pair in `ls -d 20??????_20??????`; do
  for toback in `ls $pair/*.prelb.tif 2>/dev/null`; do
    toorig=$pair/`basename $toback .prelb.tif` # e.g. geo.unw.tif
    rm $toorig
    mv $toback $toorig
  done
done
fi


# in GEOC
if [ $dogacos == 1 ]; then
# strict means gacos will be used only if available for ALL data
 #using GACOS only if there is at least for half of the epochs
 numf=`ls -d [1,2]*[0-9] | cut -d '_' -f2 | sort -u| wc -l`
 let half=$numf/2
if [ $strict == 0 ]; then
 if [ `ls ../GACOS | wc -l` -lt $half ]; then rm -r ../GACOS; dogacos=0; else dogacos=1; fi
else
 if [ `ls ../GACOS | wc -l` -lt $numf ]; then rm -r ../GACOS; dogacos=0; else dogacos=1; fi
fi
fi

# in GEOC
# backup the orig ext (and link them)
if [ $prelb_backup -gt 0 ]; then
  echo "Creating temporary backup of the input ifgs"
  #extofproc= ...
  disdir=`pwd`
  for pair in `ls -d 20??????_20??????`; do
    cd $pair
    toback=$pair.geo.$extofproc.tif
    if [ ! -f $toback ]; then
      echo "error, no original "$extofproc" tif exists for pair "$pair". removing"
      cd $disdir; mkdir -p ../GEOC.removed; if [ -d ../GEOC.removed/$pair ]; then rm -rf ../GEOC.removed/$pair; fi; mv $pair ../GEOC.removed/.
    else
      if [ ! -f $toback.prelb.tif ]; then
        mv $toback $toback.prelb.tif
        ln -s $toback.prelb.tif $toback
      else
        echo "Inconsistency detected - inform Milan for debugging"; exit
      fi
      cd $disdir
    fi
  done
fi

cd $workdir
# first store the original command:
echo $discmd > "command.in"

if [ $icams -gt 0 ]; then
  echo "regenerating ICAMS data if possible"
  if [ ! -f ~/.cdsapirc ]; then
    echo "WARNING - you did not register to cdsapi. Will now use limited access through general LiCSAR account"
    cp $LiCSAR_configpath/.cdsapirc ~/.
    #echo "ERROR - no .cdsapirc found - skipping ICAMS generation. Check with earmla (in dev)"
  fi
    create_ICAMS_frame_allepochs $frame
  #fi
  echo "Linking existing ICAMS corrections per epoch"
   mkdir -p GEOC.EPOCHS; disdir=`pwd`; cd GEOC.EPOCHS
   extfull=icams.sltd.geo.tif
   for epochpath in `ls $epochdir/20?????? -d`; do
      epoch=`basename $epochpath`
      if [ -f $epochpath/$epoch.$extfull ]; then
        if [ ! -e $epoch/$epoch.$extfull ]; then
         mkdir -p $epoch
         cd $epoch
         ln -s $epochpath/$epoch.$extfull
         cd ..
        fi
      fi
   done
fi

if [ $setides -gt 0 ]; then
  echo "checking/generating solid earth tides data"
  create_LOS_tide_frame_allepochs $frame
  disprocdir=`pwd`
  if [ $reunw -gt 0 ]; then  # in such case we correct before unwrapping
     echo "applying the SET correction"
	 # now using them to create either pha or unw tifs (to GEOC)
	 cd GEOC; disdir=`pwd`
	 hgtfile=$disdir/`ls *.geo.hgt.tif | head -n 1`
	 regt=`gmt grdinfo $hgtfile | grep registration | gawk {'print $2'}`
	 #if [ $extofproc == 'unw' ]; then grdmextra=''; else grdmextra='WRAP'; fi   # now use only for wrapped data..
	 grdmextra='WRAP'
	 for pair in `ls -d 20??????_20??????`; do
	   echo $pair
	   cd $pair
	   infile=`pwd`/$pair.geo.$extofproc.tif
	   if [ ! -L $infile ]; then
		 echo "ERROR - inconsistency detected - the file "$infile" should be already a link. Contact Milan for debugging"
		 exit
	   fi
	   date1=`echo $pair | cut -d '_' -f1`
	   date2=`echo $pair | cut -d '_' -f2`
	   outfile=`pwd`/$pair.geo.$extofproc.notides.tif
	   if [ ! -f $outfile ]; then
		 tided1=$epochdir/$date1/$date1.tide.geo.tif
		 tided2=$epochdir/$date2/$date2.tide.geo.tif   # should be A-B....
		 if [ -f $tided1 ] && [ -f $tided2 ]; then
			#echo $pair
			#
			if [ `gmt grdinfo $infile | grep registration | gawk {'print $2'}` == $regt ]; then  #Pixel ]; then
			 gmt grdmath -N $infile'=gd:Gtiff+n0' 0 NAN $tided1 $tided2 SUB 226.56 MUL SUB $grdmextra = $outfile'=gd:Gtiff'
			else
			 # half pixel issue in older frames! but ok for tides, so:
			 # echo "print('"$pair"')" >> $tmpy
			 echo "Warning, the pair "$pair" is in pixel registration. Slower workaround"
			 ifg_remove_tides.py $hgtfile $infile $tided1 $tided2 $outfile
			 # now the output is in Gridline but it says pixel (or opposite, depending on $regt)
			 # may work anyway...
			fi
			if [ -f $outfile ]; then
			  rm $infile  # only removing the link
			  ln -s `basename $outfile` `basename $infile`
			fi
		 else
		   echo "WARNING: SET estimates do not exist for pair "$pair" - perhaps one of epochs is not stored in LiCSAR_public - keeping this pair anyway"
		 fi
	   fi
	   cd $disdir
	 done
  #else   # i mean, link it anyway, as we might want to check loading to cube etc.
  else
    echo "WARNING: Without reunwrapping, the SET and iono corrs are only ready but not applied. Contact Milan - work in progress"
  fi
   # correct only on epoch level, i.e. now just link to 
   echo "Linking solid earth tide corrections per epoch"
   cd $disprocdir
   mkdir -p GEOC.EPOCHS; disdir=`pwd`; cd GEOC.EPOCHS
   extfull=tide.geo.tif
   for epochpath in `ls $epochdir/20?????? -d`; do
      epoch=`basename $epochpath`
      if [ -f $epochpath/$epoch.$extfull ]; then
        if [ ! -e $epoch/$epoch.$extfull ]; then
         mkdir -p $epoch
         cd $epoch
         ln -s $epochpath/$epoch.$extfull
         cd $disdir
        fi
      fi
   done
   cd $disprocdir
  #fi
  cd $workdir
fi


if [ $iono -gt 0 ]; then
 echo "checking/generating ionospheric correction data"
 python3 -c "from iono_correct import *; make_all_frame_epochs('"$frame"')"
 disprocdir=`pwd`
 if [ $reunw -gt 0 ]; then
	 echo "applying the ionospheric correction"
	 cd GEOC
	 # using them to either pha or unw tifs (to GEOC)
	 disdir=`pwd`
	 #hgtfile=`ls *.geo.hgt.tif | head -n 1`
	 tmpy=`pwd`/../tmp.py
	 echo "from iono_correct import correct_iono_pair;" > $tmpy
	 if [ $setides -gt 0 ]; then
		 outext=$extofproc.notides.noiono
	 else
		 outext=$extofproc.noiono
	 fi
	 for pair in `ls -d 20??????_20??????`; do
	   cd $pair
	   # here use the linked
	   infile=`pwd`/$pair.geo.$extofproc.tif
	   if [ ! -L $infile ]; then
		 echo "ERROR - inconsistency detected - the file "$infile" should be already a link. Contact Milan for debugging"
		 exit
	   fi
	   # as input, and then store as .iono.
	   # and make the link back!
	   date1=`echo $pair | cut -d '_' -f1`
	   date2=`echo $pair | cut -d '_' -f2`
	   #$epochdir
	   outfile=`pwd`/$pair.geo.$outext.tif
	   if [ ! -f $outfile ]; then
		 ionod1=$epochdir/$date1/$date1.geo.iono.code.tif
		 ionod2=$epochdir/$date2/$date2.geo.iono.code.tif   # should be A-B....
		 if [ -f $ionod1 ] && [ -f $ionod2 ]; then
			#echo $pair
			#python3 -c "from iono_correct import *;
			echo "print('"$pair"')" >> $tmpy
			echo "try:" >> $tmpy
			echo "    correct_iono_pair(frame = '"$frame"', pair = '"$pair"', ifgtype = '"$extofproc"', infile = '"$infile"', source = 'code', fixed_f2_height_km = 450, outif='"$outfile"')" >> $tmpy
			echo "except:" >> $tmpy
			echo "    print('error correcting pair "$pair"')" >> $tmpy
			#if [ $extofproc == 'unw' ]; then grdmextra=''; else grdmextra='WRAP'; fi
			#gmt grdmath $infile'=gd:Gtiff+n0' 0 NAN $ionod1 $ionod2 SUB SUB $grdmextra = $outfile'=gd:Gtiff'
			#if [ ! -f $outfile ]; then
			#  if [ -f $hgtfile ]; then
			#     echo "some error, trying to correct"
			#     gdalwarp2match.py $ionod1 $hgtfile $ionod1.tmp.tif; rm $ionod1; gdal_translate -of GTiff -co COMPRESS=DEFLATE -co PREDICTOR=3 $ionod1.tmp.tif $ionod1
			#     gdalwarp2match.py $ionod2 $hgtfile $ionod2.tmp.tif; rm $ionod2; gdal_translate -of GTiff -co COMPRESS=DEFLATE -co PREDICTOR=3 $ionod2.tmp.tif $ionod2
			#     gmt grdmath $infile'=gd:Gtiff+n0' 0 NAN $ionod1 $ionod2 SUB SUB $grdmextra = $outfile'=gd:Gtiff'
			#  fi
			#
			#rm $infile.backup 2>/dev/null
			#mv $infile $infile.backup
			#ln -s `basename $outfile` `basename $infile`
		 else
		   echo "WARNING: iono estimates do not exist for pair "$pair" - perhaps one of epochs is not stored in LiCSAR_public - keeping this pair anyway"
		 fi
	   fi
	   #if [ -f $outfile ]; then
	   #  rm $infile # that's just a link
	   #  ln -s $outfile $infile
	   #fi
	   cd $disdir
	 done
	 pairstoproc=`grep frame $tmpy | wc -l`
	 if [ $pairstoproc -gt 0 ]; then
	  echo "Correcting the ionosphere for "`grep frame $tmpy | wc -l`" pairs"
	  python3 $tmpy
	 fi
	 disdir=`pwd`
	 for pair in `ls -d 20??????_20??????`; do
	   cd $pair
	   outfile=$pair.geo.$outext.tif
	   if [ -e ${outfile} ]; then
		 # link this one instead of this link
		 ifglink=$pair.geo.$extofproc.tif
		 if [ -L $ifglink ]; then
			rm $ifglink
			ln -s $outfile $ifglink
		 else
			echo "ERROR, the file "$ifglink" should be a link - not continuing"
			exit
		 fi
	   fi
	   cd $disdir
	 done
	 rm $tmpy
  else
   echo "WARNING: Without reunwrapping, the SET and iono corrs are only ready but not applied. Contact Milan - work in progress"
  fi
  #else
   # correct only on epoch level, i.e. now just link to 
   echo "Linking iono corrections per epoch"
   cd $disprocdir
   mkdir -p GEOC.EPOCHS; disdir=`pwd`; cd GEOC.EPOCHS
   extfull=geo.iono.code.tif
   for epochpath in `ls $epochdir/20?????? -d`; do
      epoch=`basename $epochpath`
      if [ -f $epochpath/$epoch.$extfull ]; then
        if [ ! -e $epoch/$epoch.$extfull ]; then
         mkdir -p $epoch
         cd $epoch
         ln -s $epochpath/$epoch.$extfull
         cd $disdir
        fi
      fi
   done
   cd $disprocdir
  #fi
  cd $workdir
fi


#hgtfile=/gws/nopw/j04/nceo_geohazards_vol1/public/LiCSAR_products/12/012A_05443_131313/metadata/012A_05443_131313.geo.hgt.tif
 #epath=/gws/nopw/j04/nceo_geohazards_vol1/public/LiCSAR_products/12/012A_05443_131313/epochs
 #ifgspath=/gws/nopw/j04/nceo_geohazards_vol1/public/LiCSAR_products/12/012A_05443_131313/interferograms
 #frame=012A_05443_131313
 #ext=diff_pha
 #for ifg in `ls interferograms`; do
 # date1=`echo $ifg | cut -d '_' -f1`
 # date2=`echo $ifg | cut -d '_' -f2`
 # #infile=$ifgpath/$ifg/$ifg.geo.diff_pha.tif
 # infilee=$ifgspath/$ifg/$ifg.geo.diff_pha.notides.tif
 # outfilee=$ifgspath/$ifg/$ifg.geo.diff_pha.notides.noiono.tif
 # if [ ! -f $outfilee ]; then
 # if [ -f $epath/$date1/$date1.geo.iono.code.tif ] && [ -f $epath/$date2/$date2.geo.iono.code.tif ]; then
 #   #doit=1
 #   ionod1=$epath/$date1/$date1.geo.iono.code.tif
 #   ionod2=$epath/$date2/$date2.geo.iono.code.tif   # should be A-B....
 #   correct_ifg_tides_public $frame $ifg $ext
 #   gdalwarp2match.py $infilee $hgtfile $infilee'.tmp.tif'
 #   mv $infilee'.tmp.tif' $infilee
 #   gmt grdmath $infilee'=gd:Gtiff+n0' 0 NAN $ionod1 $ionod2 SUB SUB WRAP = $outfilee'=gd:Gtiff'
 #   #gmt grdmath -N $infile'=gd:Gtiff+n0' 0 NAN $tided2 $tided1 SUB -226.56 MUL SUB WRAP = $outfile'=gd:Gtiff'


if [ $reunw -gt 0 ]; then
 #echo "preparing for custom multilooking - just run ./multirun.sh"
 if [ $dogacos == 1 ]; then
  mlgeocdir=GEOCml$multi'GACOS'
 else
  mlgeocdir=GEOCml$multi
 fi
 if [ $clip == 1 ]; then
  mlgeocdir=$mlgeocdir'clip'
 fi
 mkdir $mlgeocdir
 echo cd `pwd`/$mlgeocdir > multirun.sh
 if [ $hgts == 1 ]; then
  extraparam=", hgtcorr = True"
 else
  extraparam=", hgtcorr = False"
 fi
 if [ $smooth == 1 ]; then
   echo "smooth run selected - disabling goldstein filter." # better you disable smooth, it is considered obsolete now"
  extraparam=$extraparam", smooth = True, goldstein = False"
 fi
 if [ $gammadf == 1 ]; then
   echo "using GAMMA ADF for the spatial filtering (and consistence)"
   if [ $smooth == 1 ]; then echo "WARNING: smooth option is ON - ADF will not be used. Remove -s to turn it ON."; fi
   extraparam=$extraparam", use_gamma = True"
 fi
 if [ $lowpass == 1 ]; then
  extraparam=$extraparam", lowpass = True"
 fi
 if [ $filtifg == 1 ]; then
   extraparam=$extraparam", prefer_unfiltered = False"
 fi
 # adding possibility to change consistence threshold (mask before unw) here
 extraparam=$extraparam", thres = "$thres
 if [ $keep_coh_debug == 1 ]; then
  extraparam=$extraparam", keep_coh_debug = True";
 else
  extraparam=$extraparam", keep_coh_debug = False"
 fi
 if [ $ampcoh == 1 ]; then
  extraparam=$extraparam", use_amp_coh = True";
 else
  extraparam=$extraparam", use_amp_coh = False"
 fi
 if [ $ampstab == 1 ]; then
  extraparam=$extraparam", use_amp_stab = True";
  if [ $smooth == 1 ]; then
    extraparam=$extraparam", fillby = 'gauss'";  # for test, might really help
  fi
 else
  extraparam=$extraparam", use_amp_stab = False"
 fi
 if [ $cascade == 1 ]; then
  extraparam=$extraparam", cascade = True"
  #if [ $multi -gt 2 ]; then
  # echo "setting cascade start from ML20 (10 5 3 1)"
  # extraparam=$extraparam", "
  #fi
 fi
 if [ $use_coh_stab == 1 ]; then
  extraparam=$extraparam", use_coh_stab = True"
 fi
 if [ $clip == 1 ]; then
  extraparam=$extraparam", cliparea_geo = '"$aoi"'" 
 fi
 if [ $specmag == 1 ]; then
  extraparam=$extraparam", specmag = True"
 fi
 if [ $dogacos == 1 ]; then
  extraparam=$extraparam", subtract_gacos = True" 
 fi
 if [ $defdates == 0 ]; then
  ls GEOC | grep ^20 | grep '_' | grep [0-9]\$ > pairset.txt
  extraparam=$extraparam", pairsetfile = '../pairset.txt'"
 fi
 if [ $dolocal == 1 ]; then
#  echo "your dir contains GEOC folder. using local data - not from LiCSAR_public"
  extraparam=$extraparam",  dolocal = True"
  echo "ln -s ../GEOC" >> multirun.sh
 fi
 if [ $nproc -gt 1 ]; then
   extraparam=$extraparam",  nproc = "$nproc
 fi
 if [ $landmask -lt 1 ]; then
   extraparam=$extraparam", do_landmask = False"
 fi
 if [ $outifs - gt 1 ]; then
   extraparam=$extraparam", export_to_tif = True"
 fi
 #cp $LiCSAR_procdir/$track/$frame/geo/EQA.dem_par GEOC/.
 echo "python3 -c \"from lics_unwrap import process_frame; process_frame('"$frame"', ml="$multi $extraparam")\"" >> multirun.sh
 # this seems not needed but in case of cropping, licsbas would try regenerate all missing data. so keeping this solution - may not be best if starting in local dir
 #echo "cd ..; for x in \`cat pairset.txt\`; do rm GEOC/\$x 2>/dev/null; done" >> multirun.sh
 echo "if [ \`ls [1,2]*[0-9] -d 2>/dev/null | wc -l \` -lt 2 ]; then cd ..; echo 'error processing, see processing_jasmin.* files'; exit; fi" >> multirun.sh
 echo "cd ..; cp GEOC/baselines $mlgeocdir/." >> multirun.sh
 #echo "python3 extra.py" >> multirun.sh
 #echo "python3 -c \"from LiCSAR_lib.unwrp_multiscale import process_frame; process_frame('"$frame"', ml="$multi")\"" >> multirun.sh
 #echo "cd .." >> multirun.sh
 chmod 777 multirun.sh
fi


#preparing batch file
module load $LB_version
rm -f batch_LiCSBAS.sh 2>/dev/null
copy_batch_LiCSBAS.sh >/dev/null

if [ $reunw -gt 0 ]; then # && [ $clip == 1 ]; then
 # in this case, the whole dataset should be ready for time series processing
 sed -i 's/start_step=\"01\"/start_step=\"11\"/' batch_LiCSBAS.sh
 sed -i 's/GEOCmldir=\"GEOCml${nlook}/GEOCmldir=\"'$mlgeocdir'/' batch_LiCSBAS.sh
else
 sed -i 's/start_step=\"01\"/start_step=\"02\"/' batch_LiCSBAS.sh
fi


sed -i 's/n_para=\"\"/n_para=\"'$nproc'\"/' batch_LiCSBAS.sh
sed -i 's/nlook=\"1\"/nlook=\"'$multi'\"/' batch_LiCSBAS.sh
# fix memory values
sed -i 's/p13_mem_size=\"\"/p13_mem_size=\"'$memm'\"/' batch_LiCSBAS.sh
sed -i 's/p14_mem_size=\"\"/p14_mem_size=\"'$memm'\"/' batch_LiCSBAS.sh

if [ $reunw -gt 0 ]; then
 # avoiding coh and unw cov. checking
 sed -i 's/p11_coh_thre=\"/p11_coh_thre=\"0/' batch_LiCSBAS.sh
 sed -i 's/p11_unw_thre=\"/p11_unw_thre=\"0.05/' batch_LiCSBAS.sh
 sed -i 's/p12_loop_thre=\"/p12_loop_thre=\"10/' batch_LiCSBAS.sh
 #sed -i 's/p12_loop_thre=\"/p12_loop_thre=\"30/' batch_LiCSBAS.sh   # because we would use --nullify here...
 #sed -i 's/p15_n_loop_err_thre=\"/p15_n_loop_err_thre=\"'$half'/' batch_LiCSBAS.sh
 sed -i 's/p15_resid_rms_thre=\"/p15_resid_rms_thre=\"12/' batch_LiCSBAS.sh
 sed -i 's/p15_n_loop_err_thre=\"/p15_n_loop_err_thre=\"50/' batch_LiCSBAS.sh
 sed -i 's/p15_stc_thre=\"/p15_stc_thre=\"10/' batch_LiCSBAS.sh
 #sed -i 's/start_step=\"02\"/start_step=\"16\"/' $x/batch_LiCSBAS.sh
else
 sed -i 's/p11_coh_thre=\"/p11_coh_thre=\"0.025/' batch_LiCSBAS.sh
 sed -i 's/p12_loop_thre=\"/p12_loop_thre=\"10/' batch_LiCSBAS.sh
 sed -i 's/p15_resid_rms_thre=\"/p15_resid_rms_thre=\"10/' batch_LiCSBAS.sh
# sed -i 's/p15_n_ifg_noloop_thre=\"/p15_n_ifg_noloop_thre=\"300/' batch_LiCSBAS.sh
 #sed -i 's/p15_n_loop_err_thre=\"/p15_n_loop_err_thre=\"20/' batch_LiCSBAS.sh
fi

if [ $wls == 1 ]; then
 sed -i 's/p13_inv_alg=\"\"/p13_inv_alg=\"WLS\"/' batch_LiCSBAS.sh
fi

if [ $deramp == 1 ]; then
 sed -i 's/p16_deg_deramp=\"\"/p16_deg_deramp=\"1\"/' batch_LiCSBAS.sh
fi

# set comet dev functions...
sed -i "s/^cometdev=.*/cometdev=\'"$cometdev"\'/" batch_LiCSBAS.sh


# setting those values 'everywhere' (originally it was in the modified approach):
#sed -i 's/p15_n_ifg_noloop_thre=\"/p15_n_ifg_noloop_thre=\"'$half'/' batch_LiCSBAS.sh
sed -i 's/p15_n_ifg_noloop_thre=\"/p15_n_ifg_noloop_thre=\"1000/' batch_LiCSBAS.sh   # skipping the loop closure for now - for all approaches
#sed -i 's/p16_deg_deramp=\"/p16_deg_deramp=\"1/' batch_LiCSBAS.sh   # nah, let's not deramp by default

if [ $hgtcorrlicsbas -gt 0 ]; then
 sed -i 's/p16_hgt_linear=\"n\"/p16_hgt_linear=\"y\"/' batch_LiCSBAS.sh
fi

if [ $dogacos -gt 0 ]; then
 sed -i 's/do03op_GACOS=\"n\"/do03op_GACOS=\"y\"/' batch_LiCSBAS.sh
fi

if [ $clip -gt 0 ]; then
 # lon1/lon2/lat1/lat2
 sed -i 's/do05op_clip=\"n\"/do05op_clip=\"y\"/' batch_LiCSBAS.sh
 echo $aoi | sed 's/\//\\\//g' > tmp.sedaoi
 sed -i 's/^p05_clip_range_geo=\"/p05_clip_range_geo=\"'`cat tmp.sedaoi`'/' batch_LiCSBAS.sh
 rm tmp.sedaoi
fi

if [ $maskbias -gt 0 ]; then
  sed -i 's/p15_avg_phasebias=\"\"/p15_avg_phasebias=\"1\"/' batch_LiCSBAS.sh
fi

if [ $run_jasmin -eq 1 ]; then
 echo "sending as job to JASMIN"
 rm jasmin_run.sh 2>/dev/null
 if [ $reunw -gt 0 ]; then
  cat multirun.sh > jasmin_run.sh
 fi
  #just a little export fix
  #multi=1
 #fi
 echo "module load "$LB_version >> jasmin_run.sh
 echo "./batch_LiCSBAS.sh" >> jasmin_run.sh
 
 if [ $clip -eq 1 ]; then clstr='clip'; else clstr=''; fi
 if [ $dogacos -eq 1 ]; then geocd='GEOCml'$multi"GACOS"$clstr; else geocd='GEOCml'$multi$clstr; fi
 tsdir=TS_$geocd
 if [ $reunw -eq 0 ]; then
  # so here we have already unwrapped data and we will just post-correct the ramps
  if [ $setides -gt 0 ]; then
    #echo "insert code to post-correct SET here"
    echo "Note: SET corrections will be applied to cum_filt only"
    echo "python3 -c \"from lics_tstools import *; correct_cum_from_tifs('"$tsdir"/cum_filt.h5', 'GEOC.EPOCHS', 'tide.geo.tif', 1000)\"" >> jasmin_run.sh
    #correct_cum_from_tifs(cumhdfile, tifdir = 'GEOC.EPOCHS', ext='geo.iono.code.tif', tif_scale2mm = 1, outputhdf = None, directcorrect = True)
  fi
  if [ $iono -gt 0 ]; then
    #echo "insert code to post-correct iono here"
    echo "Note: iono corrections will be applied to cum_filt only"
    echo "python3 -c \"from lics_tstools import *; correct_cum_from_tifs('"$tsdir"/cum_filt.h5', 'GEOC.EPOCHS', 'geo.iono.code.tif', 55.465/(4*np.pi))\"" >> jasmin_run.sh
  fi
 fi
 
 if [ $storeext2cube -gt 0 ]; then
  #include generation of outputs
  if [ $setides -gt 0 ]; then
    echo "Additionally, the corrections will be stored in cum.h5 as layer tide"
    echo "python3 -c \"from lics_tstools import *; correct_cum_from_tifs('"$tsdir"/cum.h5', 'GEOC.EPOCHS', 'tide.geo.tif', 1000, directcorrect = False)\"" >> jasmin_run.sh
  fi
  if [ $iono -gt 0 ]; then
    echo "Additionally, the corrections will be stored in cum.h5 as layer iono"
    echo "python3 -c \"from lics_tstools import *; correct_cum_from_tifs('"$tsdir"/cum.h5', 'GEOC.EPOCHS', 'geo.iono.code.tif', 55.465/(4*np.pi), directcorrect = False)\"" >> jasmin_run.sh
  fi
  if [ $gacos -gt 0 ]; then
    echo "Additionally, the corrections will be stored in cum.h5 as layer gacos"
    echo "python3 -c \"from lics_tstools import *; correct_cum_from_tifs('"$tsdir"/cum.h5', 'GACOS', 'sltd.geo.tif', -55.465/(4*np.pi), directcorrect = False)\"" >> jasmin_run.sh
  fi
  if [ $icams -gt 0 ]; then
    echo "Additionally, the corrections will be stored in cum.h5 as layer icams"
    echo "python3 -c \"from lics_tstools import *; correct_cum_from_tifs('"$tsdir"/cum.h5', 'GEOC.EPOCHS', 'icams.sltd.geo.tif', -55.465/(4*np.pi), directcorrect = False)\"" >> jasmin_run.sh
  fi
 fi
 echo "LiCSBAS_flt2geotiff.py -i TS_"$geocd"/results/vel.filt.mskd -p "$geocd"/EQA.dem_par -o "$frame".vel_filt.mskd.geo.tif" >> jasmin_run.sh
 echo "LiCSBAS_flt2geotiff.py -i TS_"$geocd"/results/vel.filt -p "$geocd"/EQA.dem_par -o "$frame".vel_filt.geo.tif" >> jasmin_run.sh
 echo "LiCSBAS_flt2geotiff.py -i TS_"$geocd"/results/vel.mskd -p "$geocd"/EQA.dem_par -o "$frame".vel.mskd.geo.tif" >> jasmin_run.sh
 echo "LiCSBAS_flt2geotiff.py -i TS_"$geocd"/results/vel -p "$geocd"/EQA.dem_par -o "$frame".vel.geo.tif" >> jasmin_run.sh
 echo "LiCSBAS_flt2geotiff.py -i TS_"$geocd"/results/vstd -p "$geocd"/EQA.dem_par -o "$frame".vstd.geo.tif" >> jasmin_run.sh
 echo "cp TS_"$geocd"/network/network13.png $frame'_network.png'" >> jasmin_run.sh
 echo "cp TS_"$geocd"/mask_ts.png $frame'_mask_ts.png'" >> jasmin_run.sh
 #echo "LiCSBAS_out2nc.py -i TS_GEOCml"$multi"*/cum_filt.h5 -o "$frame".nc" >> jasmin_run.sh
 echo "LiCSBAS_out2nc.py -i TS_"$geocd"/cum.h5 -o "$frame".nc" >> jasmin_run.sh
 echo "LiCSBAS_disp_img.py -i TS_"$geocd"/results/vel.filt.mskd -p "$geocd"/EQA.dem_par -c SCM.roma_r --cmin -20 --cmax 20 --kmz "$frame".vel.mskd.kmz" >> jasmin_run.sh
 echo "LiCSBAS_disp_img.py -i TS_"$geocd"/results/vel.filt.mskd -p "$geocd"/EQA.dem_par -c SCM.roma_r --cmin -20 --cmax 20 --title "$frame"_vel_filt_mskd --png "$frame".vel.mskd.png" >> jasmin_run.sh


 # jasmin proc
 cmd="bsub2slurm.sh -o processing_jasmin.out -e processing_jasmin.err -J LB_"$frame" -n "$nproc" -W "$hours":59 -M "$memmfull" -q "$que" ./jasmin_run.sh"
 echo $cmd > jasmin_run_cmd.sh
 chmod 777 jasmin_run.sh
 chmod 777 jasmin_run_cmd.sh
 #echo $cmd
 ./jasmin_run_cmd.sh
else
 nano batch_LiCSBAS.sh
 echo "now run ./batch_LiCSBAS.sh";
fi

cd $thisdir
