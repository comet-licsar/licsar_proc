#!/bin/bash
# parameters: frame [startdate] [enddate]
# e.g. 155D_02611_050400 20141001 20200205

if [ -z $1 ]; then
 echo "parameters: frame [startdate] [enddate]"
 echo "e.g. 155D_02611_050400 20141001 20200205"
 echo "parameters:"
 echo "-M 10 .... this will do extra (Gaussian-improved) multilooking and reunwrapping"
 exit
fi

multi=0
run_jasmin=1
#dogacos=0

while getopts ":M" option; do
 case "${option}" in
  M) multi=$2;
     shift
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

startdate=20141001
enddate=`date +%Y%m%d`
if [ ! -z $2 ]; then
 startdate=$2
fi
if [ ! -z $3 ]; then
 enddate=$3
fi

track=`echo $frame | cut -c -3 | sed 's/^0//' | sed 's/^0//'`

indir=$LiCSAR_public/$track/$frame/interferograms
epochdir=$LiCSAR_public/$track/$frame/epochs
metadir=$LiCSAR_public/$track/$frame/metadata
workdir=`pwd`

mkdir GEOC GACOS
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

if [ -z $2 ]; then
    #cp $epochdir/$epoch/$epoch.sltd.geo.tif ../GACOS/. 2>/dev/null
  for ifg in `ls $indir/20* -d 2>/dev/null`; do
    if [ `basename $ifg | cut -d '_' -f2` -le $enddate ]; then
      ln -s $ifg;
     fi
  done
else
 echo "nah, not ready yet, do full proc, without startdate enddate please.."
fi

#using GACOS only if there is at least for half of the files
numf=`ls | cut -d '_' -f1 | sort -u | wc -l`
let half=$numf/2
if [ `ls ../GACOS | wc -l` -lt $half ]; then rm -r ../GACOS; dogacos=0; else dogacos=1; fi

cd $workdir
if [ $multi -gt 0 ]; then
 echo "preparing for custom multilooking - just run ./multirun.sh"
 mkdir GEOCml$multi
 echo cd `pwd`/GEOCml$multi > multirun.sh
 echo "python3 -c \"from LiCSAR_lib.unwrp_multiscale import process_frame; process_frame('"$frame"', ml="$multi")\"" >> multirun.sh
 echo "cd .." >> multirun.sh
 chmod 777 multirun.sh
fi

#preparing batch file
module load LiCSBAS
copy_batch_LiCSBAS.sh >/dev/null

sed -i 's/start_step=\"01\"/start_step=\"02\"/' batch_LiCSBAS.sh
sed -i 's/n_para=\"\"/n_para=\"1\"/' batch_LiCSBAS.sh
if [ $multi -gt 0 ]; then
 sed -i 's/nlook=\"1\"/nlook=\"'$multi'\"/' batch_LiCSBAS.sh
 sed -i 's/p12_loop_thre=\"/p12_loop_thre=\"10/' batch_LiCSBAS.sh
 sed -i 's/p15_n_ifg_noloop_thre=\"/p15_n_ifg_noloop_thre=\"300/' batch_LiCSBAS.sh
 sed -i 's/p15_n_loop_err_thre=\"/p15_n_loop_err_thre=\"200/' batch_LiCSBAS.sh
 sed -i 's/p15_resid_rms_thre=\"/p15_resid_rms_thre=\"50/' batch_LiCSBAS.sh
 #sed -i 's/start_step=\"02\"/start_step=\"16\"/' $x/batch_LiCSBAS.sh
 sed -i 's/p16_deg_deramp=\"/p16_deg_deramp=\"2/' batch_LiCSBAS.sh
fi
if [ $dogacos -gt 0 ]; then
 sed -i 's/do03op_GACOS=\"n\"/do03op_GACOS=\"y\"/' batch_LiCSBAS.sh
fi

if [ $run_jasmin -eq 1 ]; then
 echo "sending as job to JASMIN"
 rm jasmin_run.sh 2>/dev/null
 if [ $multi -gt 0 ]; then
  cat multirun.sh > jasmin_run.sh
 fi
 echo "module load LiCSBAS" >> jasmin_run.sh
 echo "./batch_LiCSBAS.sh" >> jasmin_run.sh
 #include generation of outputs
 echo "LiCSBAS_flt2geotiff.py -i TS_GEOCml"$multi"/results/vel.filt.mskd -p GEOCml"$multi"/EQA.dem_par -o "$frame".vel_deramp.mskd.geo.tif" >> jasmin_run.sh
 echo "LiCSBAS_flt2geotiff.py -i TS_GEOCml"$multi"/results/vel.filt -p GEOCml"$multi"/EQA.dem_par -o "$frame".vel_deramp.geo.tif" >> jasmin_run.sh
 echo "LiCSBAS_flt2geotiff.py -i TS_GEOCml"$multi"/results/vel.mskd -p GEOCml"$multi"/EQA.dem_par -o "$frame".vel.mskd.geo.tif" >> jasmin_run.sh
 echo "LiCSBAS_flt2geotiff.py -i TS_GEOCml"$multi"/results/vel -p GEOCml"$multi"/EQA.dem_par -o "$frame".vel.geo.tif" >> jasmin_run.sh
 echo "LiCSBAS_flt2geotiff.py -i TS_GEOCml"$multi"/results/vstd -p GEOCml"$multi"/EQA.dem_par -o "$frame".vstd.geo.tif" >> jasmin_run.sh
 cmd="bsub2slurm.sh -o processing_jasmin.out -e processing_jasmin.err -J LB_"$frame" -n 1 -W 14:00 -q comet ./jasmin_run.sh"
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
