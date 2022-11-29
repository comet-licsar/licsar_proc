#!/bin/bash

echo "-----------------------------------------"
echo "Welcome to the JASMIN residence of LiCSAR"
echo "-----------------------------------------"
echo ""

sleep 1

echo "setting new default path to BATCH_CACHE_DIR and creating it: "
source ~/.bashrc
sed -i 's/export BATCH_CACHE_DIR/\#export BATCH_CACHE_DIR/' ~/.bashrc



if [ `grep -c "export oldbatch=" ~/.bashrc` -lt 1 ]; then
 oldbatch=$BATCH_CACHE_DIR
 if [ ! -z $oldbatch ]; then
  echo "export oldbatch="$oldbatch >> ~/.bashrc
 fi
 echo "export BATCH_CACHE_DIR=/work/scratch-pw3/licsar/"$USER"/batchdir" >> ~/.bashrc
fi


export BATCH_CACHE_DIR=/work/scratch-pw3/licsar/$USER/batchdir
export LiCSAR_temp=/work/scratch-pw3/licsar/$USER/LiCSAR_temp

echo "\$BATCH_CACHE_DIR="$BATCH_CACHE_DIR
echo "\$LiCSAR_temp="$LiCSAR_temp

mkdir -p $BATCH_CACHE_DIR
mkdir -p $LiCSAR_temp

if [ -z $oldbatch ]; then
 echo "no previous BATCH_CACHE_DIR was set (this is ok for New Users)"
else
if [ $oldbatch == $BATCH_CACHE_DIR ]; then
 echo "ERROR - it seems you might have run this procedure already - cancelling"
 exit
 #"if so, please run following command manually:"
 #echo "rsync -r YOUROLDBATCHDIR/* "$BATCH_CACHE_DIR
else
 cd $oldbatch;
 echo "and now the long part - copying your current batchdir data to the new disk"
 for x in `ls`; do
  echo "copying "$x;
  rsync -r $x $BATCH_CACHE_DIR 2>/dev/null
 done
fi
fi

echo "umask 002" >> ~/.bashrc
echo "umask 002" >> ~/.profile

echo "all done - thank you"
echo "(you may check the \$BATCH_CACHE_DIR where all your processing will be performed)"
echo "REMEMBER data stored in this folder (on scratch disk) are auto-deleted after 28 days"
echo "so your task is to process data, inform us once done so we can store them, copy only necessary data to your work directory in vol1 or vol2 disk"
