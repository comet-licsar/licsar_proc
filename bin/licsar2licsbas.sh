#!/bin/bash
# parameters: frame [startdate] [enddate]
# e.g. 155D_02611_050400 20141001 20200205

if [ -z $1 ]; then
 echo "parameters: frame [startdate] [enddate]"
 echo "e.g. 155D_02611_050400 20141001 20200205"
 exit
fi

frame=$1

startdate=20141001
enddate=`date +%Y%m%d`
if [ ! -z $2 ]; then
 startdate=$2
fi
if [ ! -z $3 ]; then
 enddate=$3
fi

track=`echo $frame | cut -c -3 | sed 's/^0//' | sed 's/^0//'`

indir=$LiCSAR_public/$track/$frame/products
metadir=$LiCSAR_public/$track/$frame/metadata
workdir=`pwd`

mkdir GEOC
cd GEOC
for meta in E N U hgt; do
 ln -s $metadir/$frame.geo.$meta.tif
done
ln -s $metadir/baselines

source $metadir/metadata.txt #this will bring 'master=' info
ln -s $indir/$master/$master.geo.mli.tif

for epoch in `ls $indir/epochs`; do
  if [ $epoch -ge $startdate ] && [ $epoch -lt $enddate ]; then
    for ifg in `ls $indir/$epoch* -d 2>/dev/null`; do
     if [ `basename $ifg | cut -d '_' -f2` -le $enddate ]; then
      ln -s $ifg;
     fi
    done
  fi
done

cd $workdir

#preparing batch file
module load LiCSBAS
copy_batch_LiCSBAS.sh

sed -i 's/start_step=\"01\"/start_step=\"02\"/' batch_LiCSBAS.sh
nano batch_LiCSBAS.sh
./batch_LiCSBAS.sh
