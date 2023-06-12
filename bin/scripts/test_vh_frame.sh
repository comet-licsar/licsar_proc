#!/bin/bash

echo "one off test for VH data"
frame=050D_05049_060600
sdate=20141001
#sdate=20230101
edate=20230520

#mkdir -p $BATCH_CACHE_DIR/VHTEST
#cd $BATCH_CACHE_DIR/VHTEST
#origba=$BATCH_CACHE_DIR
#export BATCH_CACHE_DIR=`pwd`

if [ -d $BATCH_CACHE_DIR/$frame ]; then
echo "careful, this frame already exists"
mv $BATCH_CACHE_DIR/$frame/RSLC $BATCH_CACHE_DIR/$frame/RSLC.orig
mv $BATCH_CACHE_DIR/$frame/SLC $BATCH_CACHE_DIR/$frame/SLC.orig
fi

licsar_make_frame.sh -n -f -P $frame 1 1 ${sdate:0:4}-${sdate:4:2}-${sdate:6:2} ${edate:0:4}-${edate:4:2}-${edate:6:2}
touch $BATCH_CACHE_DIR/$frame/lmf_locked
touch $BATCH_CACHE_DIR/$frame/lmf_locked

cd $BATCH_CACHE_DIR/$frame
# remove the master date
m=`get_master`
rm -rf RSLC/* SLC/$m

cd ..; mkdir VHTEST
mv $frame VHTEST/.
cd VHTEST/$frame
rm -rf GEOC GEOC.MLI IFG

# generate all VH SLCs, make sure master date is included!
LiCSAR_01_mk_images.py -n -f $frame -d . -s $m -e $m -a 4 -r 20 -x
if [ ! -d SLC/$m ]; then
 echo "ERROR - you need to download zips for the reference data first"
 exit
fi
# make master rslcs
mkdir RSLC/$m
for x in `ls SLC/$m`; do
 ln -s `pwd`/SLC/$m/$x `pwd`/RSLC/$m/`echo $x | sed 's/slc/rslc/'`
done

# extract all requested SLCs that have LUTs
ls LUT > luts.txt
LiCSAR_01_mk_images.py -n -f $frame -d . -l luts.txt -s $sdate -e $edate -a 4 -r 20 -x

# run coreg (will use recoreg with LUTs)
LiCSAR_02_coreg.py -f $frame -d . -i -m $m

echo "done"
pwd
