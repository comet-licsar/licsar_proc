#!/bin/bash

echo "this script might take long time to run - please run in tmux, i.e. do:"
echo "module load licsar_framebatch"
echo "source /gws/smf/j04/nceo_geohazards/software/mambalics/load_mambalics.rc"
echo "tmux"
echo ""
echo "if you are already in tmux, just go on - if not, cancel by pressing CTRL-C"
echo "(waiting 10 seconds before really starting)"

sleep 10
echo "ok, starting the procedure"

sleep 2

echo "setting new default path to BATCH_CACHE_DIR and creating it: "
source ~/.bashrc
oldbatch=$BATCH_CACHE_DIR
sed -i 's/export BATCH_CACHE_DIR/\#export BATCH_CACHE_DIR/' ~/.bashrc
echo "oldbatch="$oldbatch >> ~/.bashrc

export BATCH_CACHE_DIR=/work/scratch-pw3/licsar/$USER/batchdir
export LiCSAR_temp=/work/scratch-pw3/licsar/$USER/LiCSAR_temp

echo $BATCH_CACHE_DIR
echo $LiCSAR_temp

mkdir -p $BATCH_CACHE_DIR
mkdir -p $LiCSAR_temp

echo "and now the long part - copying your current batchdir data to the new disk"
if [ $oldbatch == $BATCH_CACHE_DIR ]; then
 echo "ERROR - it seems you might have run this procedure already - if so, please run following command manually:"
 echo "rsync -r YOUROLDBATCHDIR/* "$BATCH_CACHE_DIR
else
 cd $oldbatch;
 for x in `ls`; do
  echo "copying "$x;
  rsync -r $x $BATCH_CACHE_DIR
 done
fi

echo "all done - thank you"
