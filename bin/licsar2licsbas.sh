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
 echo "-s ....... use coherence stability index instead of orig coh per ifg (experimental - might help against loop closure errors, maybe)"
 echo "-k ....... use cohratio everywhere (i.e. for unwrapping, rather than orig coh - this is experimental attempt)"
 echo "-H ....... this will use hgt to support unwrapping (only if using reunwrapping)"
 echo "-T ....... use testing version of LiCSBAS"
 echo "-S ....... strict mode - e.g. in case of GACOS, use it only if available for ALL ifgs"
 echo "-G lon1/lon2/lat1/lat2  .... clip to this AOI"
 echo "----------------"
 #echo "note: in case you combine -G and -u, the result will be in clip folder without GACOS! (still not smoothly combined reunw->licsbas, todo!)"  # updated on 2022-04-07
 #echo "(note: if you do -M 1, it will go for reprocessing using the cascade/multiscale unwrap approach - in testing, please give feedback to Milan)"
 exit
fi

dolocal=0
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
#LB_version=licsbas_comet_dev
#LB_version=LiCSBAS_testing

while getopts ":M:HucTsSkG:" option; do
 case "${option}" in
  M) multi=${OPTARG};
     #shift
     ;;
  H) hgts=1;
     #shift
     ;;
  u) reunw=1;
     #shift
     ;;
  c) cascade=1;
     ;;
  s) use_coh_stab=1;
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
 esac
done
shift $((OPTIND -1))



frame=$1
if [ `echo $frame | grep -c '_'` -lt 1 ]; then
 echo "this is not a frame - check your input parameters please"
fi

thisdir=`pwd`
if [ ! `pwd | rev | cut -d '/' -f1 | rev` == $frame ]; then
 mkdir $frame
 cd $frame
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

if [ -d GEOC ]; then
 echo "warning - GEOC folder detected. will use its contents for processing, rather than link from LiCSAR_public"
 dolocal=1;
fi

mkdir GEOC GACOS 2>/dev/null
cd GEOC
for meta in E N U hgt; do
 ln -s $metadir/$frame.geo.$meta.tif
done
ln -s $metadir/baselines

source $metadir/metadata.txt #this will bring 'master=' info
if [ ! -f $epochdir/$master/$master.geo.mli.tif ]; then
 echo "regenerating missing master MLI geotiff"
 create_geoctiffs_to_pub.sh -M $LiCSAR_procdir/$track/$frame $master
 mkdir -p $epochdir/$master
 mv $LiCSAR_procdir/$track/$frame/GEOC.MLI/$master/* $epochdir/$master/.
 rmdir $LiCSAR_procdir/$track/$frame/GEOC.MLI/$master
fi
if [ -f $epochdir/$master/$master.geo.mli.tif ]; then
 ln -s $epochdir/$master/$master.geo.mli.tif
else
 echo "warning, primary epoch not in public dir. trying to use other, expect possible size issues.."
 lastimg=`ls $epochdir/*/*.geo.mli.tif | tail -n1`
 ln -s $lastimg
fi

for epoch in `ls $epochdir`; do
  if [ $epoch -ge $startdate ] && [ $epoch -le $enddate ]; then
    gacosfile=$epochdir/$epoch/$epoch.sltd.geo.tif
    if [ -f $gacosfile ]; then
     ln -s $gacosfile $workdir/GACOS/$epoch.sltd.geo.tif
    fi
  fi
done

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

# strict means gacos will be used only if available for ALL data
 #using GACOS only if there is at least for half of the files
 numf=`ls -d 20*[0-9] | cut -d '_' -f2 | sort -u| wc -l`
 let half=$numf/2
if [ $strict == 0 ]; then
 if [ `ls ../GACOS | wc -l` -lt $half ]; then rm -r ../GACOS; dogacos=0; else dogacos=1; fi
else
 if [ `ls ../GACOS | wc -l` -lt $numf ]; then rm -r ../GACOS; dogacos=0; else dogacos=1; fi
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
 if [ $keep_coh_debug == 1 ]; then
  extraparam=$extraparam", keep_coh_debug = True";
 else
  extraparam=$extraparam", keep_coh_debug = False"
 fi
 if [ $cascade == 1 ]; then
  extraparam=$extraparam", cascade = True"
 fi
 if [ $use_coh_stab == 1 ]; then
  extraparam=$extraparam", use_coh_stab = True"
 fi
 if [ $clip == 1 ]; then
  extraparam=$extraparam", cliparea_geo = '"$aoi"'" 
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
 #cp $LiCSAR_procdir/$track/$frame/geo/EQA.dem_par GEOC/.
 echo "python3 -c \"from LiCSAR_lib.unwrp_multiscale import process_frame; process_frame('"$frame"', ml="$multi $extraparam")\"" >> multirun.sh
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


sed -i 's/n_para=\"\"/n_para=\"1\"/' batch_LiCSBAS.sh
sed -i 's/nlook=\"1\"/nlook=\"'$multi'\"/' batch_LiCSBAS.sh


if [ $reunw -gt 0 ]; then
 sed -i 's/p12_loop_thre=\"/p12_loop_thre=\"10/' batch_LiCSBAS.sh
 #sed -i 's/p12_loop_thre=\"/p12_loop_thre=\"30/' batch_LiCSBAS.sh   # because we would use --nullify here...
 sed -i 's/p15_n_loop_err_thre=\"/p15_n_loop_err_thre=\"'$half'/' batch_LiCSBAS.sh
 sed -i 's/p15_resid_rms_thre=\"/p15_resid_rms_thre=\"5/' batch_LiCSBAS.sh
 #sed -i 's/start_step=\"02\"/start_step=\"16\"/' $x/batch_LiCSBAS.sh
else
 sed -i 's/p11_coh_thre=\"/p11_coh_thre=\"0.025/' batch_LiCSBAS.sh
 sed -i 's/p12_loop_thre=\"/p12_loop_thre=\"10/' batch_LiCSBAS.sh
 sed -i 's/p15_resid_rms_thre=\"/p15_resid_rms_thre=\"10/' batch_LiCSBAS.sh
# sed -i 's/p15_n_ifg_noloop_thre=\"/p15_n_ifg_noloop_thre=\"300/' batch_LiCSBAS.sh
 sed -i 's/p15_n_loop_err_thre=\"/p15_n_loop_err_thre=\"20/' batch_LiCSBAS.sh
fi

# setting those values 'everywhere' (originally it was in the modified approach):
sed -i 's/p15_n_ifg_noloop_thre=\"/p15_n_ifg_noloop_thre=\"'$half'/' batch_LiCSBAS.sh
sed -i 's/p16_deg_deramp=\"/p16_deg_deramp=\"1/' batch_LiCSBAS.sh
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
 echo "LiCSBAS_flt2geotiff.py -i TS_GEOCml"$multi"*/results/vel.filt.mskd -p "$geocd"/EQA.dem_par -o "$frame".vel_deramp.mskd.geo.tif" >> jasmin_run.sh
 echo "LiCSBAS_flt2geotiff.py -i TS_GEOCml"$multi"*/results/vel.filt -p "$geocd"/EQA.dem_par -o "$frame".vel_deramp.geo.tif" >> jasmin_run.sh
 echo "LiCSBAS_flt2geotiff.py -i TS_GEOCml"$multi"*/results/vel.mskd -p "$geocd"/EQA.dem_par -o "$frame".vel.mskd.geo.tif" >> jasmin_run.sh
 echo "LiCSBAS_flt2geotiff.py -i TS_GEOCml"$multi"*/results/vel -p "$geocd"/EQA.dem_par -o "$frame".vel.geo.tif" >> jasmin_run.sh
 echo "LiCSBAS_flt2geotiff.py -i TS_GEOCml"$multi"*/results/vstd -p "$geocd"/EQA.dem_par -o "$frame".vstd.geo.tif" >> jasmin_run.sh
 echo "cp TS_GEOCml"$multi"*/network/network13.png $frame'_network.png'" >> jasmin_run.sh
 echo "cp TS_GEOCml"$multi"*/mask_ts.png $frame'_mask_ts.png'" >> jasmin_run.sh
 #echo "LiCSBAS_out2nc.py -i TS_GEOCml"$multi"*/cum_filt.h5 -o "$frame".nc" >> jasmin_run.sh
 echo "LiCSBAS_out2nc.py -i TS_GEOCml"$multi"*/cum.h5 -o "$frame".nc" >> jasmin_run.sh
 hours=14
 if [ $multi -eq 1 ]; then
  hours=23
 fi
 cmd="bsub2slurm.sh -o processing_jasmin.out -e processing_jasmin.err -J LB_"$frame" -n 1 -W "$hours":00 -M 8192 -q comet ./jasmin_run.sh"
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
