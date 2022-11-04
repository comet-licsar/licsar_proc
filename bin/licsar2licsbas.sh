#!/bin/bash
# parameters: frame [startdate] [enddate]
# e.g. 155D_02611_050400 20141001 20200205

if [ -z $1 ]; then
 echo "parameters: frame [startdate] [enddate]"
 echo "e.g. 155D_02611_050400 20141001 20200205"
 echo "parameters:"
 echo "-M 10 .... this will do extra multilooking (in this example, 10x multilooking)"
 echo "-u ....... use the (extra Gaussian-improved multilooking and) reunwrapping procedure (useful if multilooking..)"
 echo "-c ....... if the reunwrapping is to be performed, use cascade (might be better, especially when with shores)"
 echo "-l ....... if the reunwrapping is to be performed, would do lowpass filter (should be safe unless in tricky areas as islands; good to use by default)"
 echo "-P ....... prioritise, i.e. use comet queue instead of short-serial"
 echo "-n 1 ..... number of processors (by default: 1, used also for reunwrapping, although not tested well yet)"
 #echo "-C ....... use coherence stability index instead of orig coh per ifg (experimental - might help against loop closure errors, maybe)"
 #echo "-k ....... use cohratio everywhere (i.e. for unwrapping, rather than orig coh - this is experimental attempt)"
 echo "-H ....... this will use hgt to support unwrapping (only if using reunwrapping)"
 echo "-m ....... with reunwrapping, use coh based on spectral magnitude (otherwise nyquist-limited phase difference coherence) - recommended param"
 echo "-T ....... use testing version of LiCSBAS"
 echo "-d ....... use the dev parameters for the testing version of LiCSBAS (currently: this will use --fast, --nopngs and --nullify)"
 echo "-t 0.35 ... change coherence threshold to 0.35 (default) during reunwrapping (-u)"
 echo "-g ....... use GACOS if available - NOTE THIS WAS ON BY DEFAULT TILL SEP 2022, BUT NOT ANYMORE"
 echo "-G lon1/lon2/lat1/lat2  .... clip to this AOI"
 echo "-W ....... use WLS for the inversion (coherence-based)"
 echo "----------------"
 echo "some older (not recommended anymore) parameters:"
 echo "-s ....... if the reunwrapping is to be performed, use smooth operation before adding back residuals (could be wrong - rather use lowpass+Goldstein)"
 echo "-S ....... strict mode - e.g. in case of GACOS, use it only if available for ALL ifgs"
 echo "(-R ....... prioritise through comet responder)"
 #echo "note: in case you combine -G and -u, the result will be in clip folder without GACOS! (still not smoothly combined reunw->licsbas, todo!)"  # updated on 2022-04-07
 #echo "(note: if you do -M 1, it will go for reprocessing using the cascade/multiscale unwrap approach - in testing, please give feedback to Milan)"
 exit
fi

thres=0.35
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
LB_version=LiCSBAS
cascade=0
smooth=0
lowpass=0
wls=0
cometdev=0
specmag=0
nproc=1
que='short-serial'
#LB_version=licsbas_comet_dev
#LB_version=LiCSBAS_testing

while getopts ":M:HucTsdSClWgmPRkG:t:n:" option; do
 case "${option}" in
  M) multi=${OPTARG};
     #shift
     ;;
  n) nproc=${OPTARG};
     #que='par-single'; # unless changed to comet queue
     ;;
  H) hgts=1;
     #shift
     ;;
  m) specmag=1;
     ;;
  u) reunw=1;
     #shift
     ;;
  c) cascade=1;
     ;;
  g) dogacos=1;
     ;;
  d) cometdev=1;
     ;;
  C) use_coh_stab=1;
     #shift
     ;;
  s) smooth=1;
     #shift
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
  T) LB_version=licsbas_comet_dev;
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
 if [ que == 'short-serial' ]; then
  que='par-single';
 fi 
fi

frame=$1
if [ `echo $frame | grep -c '_'` -lt 1 ]; then
 echo "this is not a frame - check your input parameters please, yet continuing"
fi

if [ -d GEOC ]; then
 echo "warning - GEOC folder detected. will use its contents for processing, rather than link from LiCSAR_public"
 dolocal=1;
fi

thisdir=`pwd`
if [ $dolocal == 0 ]; then
if [ ! `pwd | rev | cut -d '/' -f1 | rev` == $frame ]; then
 mkdir $frame
 cd $frame
fi
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
metadir=$LiCSAR_public/$track/$frame/metadata
workdir=`pwd`
source $metadir/metadata.txt #this will bring 'master=' info


 # setting for JASMIN processing
 hours=14
 if [ $multi -eq 1 ]; then
  hours=23
 fi
 memm=8192 # requesting 8GB RAM per processor
 let memm=$nproc'*'$memm
 
 
mkdir GEOC 2>/dev/null
if [ $dogacos == 1 ]; then
  mkdir GACOS 2>/dev/null
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
  echo "warning, doing local proc only - avoiding possible problems by just removing ENU files"
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
     ln -s $gacosfile $workdir/GACOS/$epoch.sltd.geo.tif
    fi
  fi
done
fi

if [ $dolocal == 0 ]; then
if [ ! -z $2 ]; then
  echo "limiting the dataset to dates between "$startdate" and "$enddate
    #cp $epochdir/$epoch/$epoch.sltd.geo.tif ../GACOS/. 2>/dev/null
  for ifg in `ls $indir/20* -d 2>/dev/null`; do
   if [ `basename $ifg | cut -d '_' -f1` -ge $startdate ]; then
    if [ `basename $ifg | cut -d '_' -f2` -le $enddate ]; then
      ln -s $ifg;
    fi
   fi
  done
else
 for ifg in `ls $indir/20* -d 2>/dev/null`; do ln -s $ifg; done
 #echo "nah, not ready yet, do full proc, without startdate enddate please.."
fi
fi


if [ $dogacos == 1 ]; then
# strict means gacos will be used only if available for ALL data
 #using GACOS only if there is at least for half of the files
 numf=`ls -d 20*[0-9] | cut -d '_' -f2 | sort -u| wc -l`
 let half=$numf/2
if [ $strict == 0 ]; then
 if [ `ls ../GACOS | wc -l` -lt $half ]; then rm -r ../GACOS; dogacos=0; else dogacos=1; fi
else
 if [ `ls ../GACOS | wc -l` -lt $numf ]; then rm -r ../GACOS; dogacos=0; else dogacos=1; fi
fi
fi

cd $workdir
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
   echo "smooth run selected - disabling goldstein filter. better you disable smooth, it is considered obsolete now"
  extraparam=", smooth = True, goldstein = False"
 fi
 if [ lowpass == 1 ]; then
  extraparam=", lowpass = True"
 fi
 # adding possibility to change coh threshold here
 extraparam=", thres = "$thres
 if [ $keep_coh_debug == 1 ]; then
  extraparam=$extraparam", keep_coh_debug = True";
 else
  extraparam=$extraparam", keep_coh_debug = False"
 fi
 if [ $cascade == 1 ]; then
  extraparam=$extraparam", cascade = True"
  if [ $multi -gt 2 ]; then
   echo "setting cascade start from ML20 (10 5 3 1)"
   extraparam=$extraparam", "
  fi
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
  ls GEOC | grep ^20 | grep '_' > pairset.txt
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
 #cp $LiCSAR_procdir/$track/$frame/geo/EQA.dem_par GEOC/.
 echo "python3 -c \"from LiCSAR_lib.lics_unwrap import process_frame; process_frame('"$frame"', ml="$multi $extraparam")\"" >> multirun.sh
 # this seems not needed but in case of cropping, licsbas would try regenerate all missing data. so keeping this solution - may not be best if starting in local dir
 #echo "cd ..; for x in \`cat pairset.txt\`; do rm GEOC/\$x 2>/dev/null; done" >> multirun.sh
 echo "cd ..; cp GEOC/baselines $mlgeocdir/." >> multirun.sh
 #echo "python3 extra.py" >> multirun.sh
 #echo "python3 -c \"from LiCSAR_lib.unwrp_multiscale import process_frame; process_frame('"$frame"', ml="$multi")\"" >> multirun.sh
 #echo "cd .." >> multirun.sh
 chmod 777 multirun.sh
fi


#preparing batch file
module load $LB_version
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
 sed -i 's/p12_loop_thre=\"/p12_loop_thre=\"10/' batch_LiCSBAS.sh
 #sed -i 's/p12_loop_thre=\"/p12_loop_thre=\"30/' batch_LiCSBAS.sh   # because we would use --nullify here...
 #sed -i 's/p15_n_loop_err_thre=\"/p15_n_loop_err_thre=\"'$half'/' batch_LiCSBAS.sh
 sed -i 's/p15_resid_rms_thre=\"/p15_resid_rms_thre=\"12/' batch_LiCSBAS.sh
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

# set comet dev functions...
sed -i "s/^cometdev=.*/cometdev=\'"$cometdev"\'/" batch_LiCSBAS.sh


# setting those values 'everywhere' (originally it was in the modified approach):
#sed -i 's/p15_n_ifg_noloop_thre=\"/p15_n_ifg_noloop_thre=\"'$half'/' batch_LiCSBAS.sh
sed -i 's/p15_n_ifg_noloop_thre=\"/p15_n_ifg_noloop_thre=\"1000/' batch_LiCSBAS.sh   # skipping the loop closure for now - for all approaches
#sed -i 's/p16_deg_deramp=\"/p16_deg_deramp=\"1/' batch_LiCSBAS.sh   # nah, let's not deramp by default
sed -i 's/p16_hgt_linear=\"n\"/p16_hgt_linear=\"y\"/' batch_LiCSBAS.sh

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
 #include generation of outputs
 if [ $clip -eq 1 ]; then clstr='clip'; else clstr=''; fi
 if [ $dogacos -eq 1 ]; then geocd='GEOCml'$multi"GACOS"$clstr; else geocd='GEOCml'$multi$clstr; fi
 echo "LiCSBAS_flt2geotiff.py -i TS_"$geocd"/results/vel.filt.mskd -p "$geocd"/EQA.dem_par -o "$frame".vel_filt.mskd.geo.tif" >> jasmin_run.sh
 echo "LiCSBAS_flt2geotiff.py -i TS_"$geocd"/results/vel.filt -p "$geocd"/EQA.dem_par -o "$frame".vel_filt.geo.tif" >> jasmin_run.sh
 echo "LiCSBAS_flt2geotiff.py -i TS_"$geocd"/results/vel.mskd -p "$geocd"/EQA.dem_par -o "$frame".vel.mskd.geo.tif" >> jasmin_run.sh
 echo "LiCSBAS_flt2geotiff.py -i TS_"$geocd"/results/vel -p "$geocd"/EQA.dem_par -o "$frame".vel.geo.tif" >> jasmin_run.sh
 echo "LiCSBAS_flt2geotiff.py -i TS_"$geocd"/results/vstd -p "$geocd"/EQA.dem_par -o "$frame".vstd.geo.tif" >> jasmin_run.sh
 echo "cp TS_"$geocd"/network/network13.png $frame'_network.png'" >> jasmin_run.sh
 echo "cp TS_"$geocd"/mask_ts.png $frame'_mask_ts.png'" >> jasmin_run.sh
 #echo "LiCSBAS_out2nc.py -i TS_GEOCml"$multi"*/cum_filt.h5 -o "$frame".nc" >> jasmin_run.sh
 echo "LiCSBAS_out2nc.py -i TS_"$geocd"/cum_filt.h5 -o "$frame".nc" >> jasmin_run.sh


 # jasmin proc
 cmd="bsub2slurm.sh -o processing_jasmin.out -e processing_jasmin.err -J LB_"$frame" -n "$nproc" -W "$hours":00 -M "$memm" -q "$que" ./jasmin_run.sh"
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
