#!/bin/bash
source $LiCSARpath/lib/LiCSAR_bash_lib.sh

echo "one off test for VH data"
echo "it will autodownload missing data in given time range, and then process VH layers of only epochs for which we have LUT files already."

if [ -z $1 ]; then
 echo "parameters: frame [DATEIN] [DATEOUT]"
 echo "(use dates as e.g 20141001 20230512)"
 exit
fi
frame=$1
sdate=20141001
edate=20290520
if [ ! -z $2 ]; then sdate=$2; fi
if [ ! -z $3 ]; then edate=$3; fi
#mkdir -p $BATCH_CACHE_DIR/VHTEST
#cd $BATCH_CACHE_DIR/VHTEST
#origba=$BATCH_CACHE_DIR
#export BATCH_CACHE_DIR=`pwd`

if [ -d $BATCH_CACHE_DIR/$frame ]; then
echo "careful, this frame already exists. but it should work anyway fine"
#mv $BATCH_CACHE_DIR/$frame/RSLC $BATCH_CACHE_DIR/$frame/RSLC.orig
#mv $BATCH_CACHE_DIR/$frame/SLC $BATCH_CACHE_DIR/$frame/SLC.orig
fi

if [ -d $BATCH_CACHE_DIR/$frame/VH ]; then
 echo "a previous VH processing detected. might work anyway but not tested (maybe delete/rename that VH dir)"
fi
licsar_make_frame.sh -n -f -P $frame 1 1 ${sdate:0:4}-${sdate:4:2}-${sdate:6:2} ${edate:0:4}-${edate:4:2}-${edate:6:2}
#touch $BATCH_CACHE_DIR/$frame/lmf_locked

cd $BATCH_CACHE_DIR/$frame
m=`get_master`
echo "checking/extracting all LUT available tables for given dates"
mkdir -p LUT
cd LUT
tr=`track_from_frame $frame`
for lut in `ls $LiCSAR_procdir/$tr/$frame/LUT`; do
 datee=`echo $lut | cut -d '.' -f1`
 if [ $datee -gt $sdate ] && [ $datee -lt $edate ]; then
   if [ ! -d $datee ]; then
    echo "extracting LUT for "$datee
    7za x $LiCSAR_procdir/$tr/$frame/LUT/$lut >/dev/null
   fi
 fi
done
cd ..
mkdir -p VH; cd VH; ln -s ../geo; ln -s ../LUT; cp ../$frame* .; mkdir SLC RSLC tab log LOGS
# remove all *LCs, including the master date
#rm -rf RSLC/* SLC/* #$m

#cd ..; mkdir VHTEST
#mv $frame VHTEST/.
#cd VHTEST/$frame
#rm -rf GEOC GEOC.MLI IFG

# generate all VH SLCs, starting by the master date
framebatch_data_refill.sh $frame ${m:0:4}-${m:4:2}-${m:6:2} ${m:0:4}-${m:4:2}-${m:6:2} # will skip the M but at least get its filenames
cd $LiCSAR_SLC
# this is to download M images as they would probably not get through the NLA request
for todown in `cat $BATCH_CACHE_DIR/$frame/$frame'_scihub.list' | rev | cut -d '/' -f 1 | rev`; do 
 wget_alaska $todown
 arch2DB.py -f `pwd`/$todown
done
cd -
LiCSAR_01_mk_images.py -n -f $frame -d . -s $m -e $m -a 4 -r 20 -x
if [ ! -d SLC/$m ]; then
 #echo "ERROR - you need to download zips for the reference data first"
 echo "ERROR preparing VH of the reference epoch, cannot continue"
 exit
fi
# make master rslcs
#mkdir RSLC/$m
#for x in `ls SLC/$m`; do
# ln -s `pwd`/SLC/$m/$x `pwd`/RSLC/$m/`echo $x | sed 's/slc/rslc/'`
#done

# extract all requested SLCs that have LUTs
ls LUT > luts.txt
LiCSAR_01_mk_images.py -n -f $frame -d . -l luts.txt -s $sdate -e $edate -a 4 -r 20 -x > vh_mk_images.out 2> vh_mk_images.err

# run coreg (will use recoreg with LUTs)
LiCSAR_02_coreg.py -f $frame -d . -i -m $m > vh_coreg.out 2> vh_coreg.err

#rm $BATCH_CACHE_DIR/$frame/lmf_locked  #not there anymore

echo "done, see following directory: "
pwd
