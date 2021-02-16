#/bin/bash
#to be run from frame directory!

frame=`pwd | rev | cut -d '/' -f1 | rev`
if [ ! `echo $frame | cut -d '_' -f3 | wc -c` -eq 7 ]; then echo "to be run from frame directory"; exit; fi

tr=`echo $frame | cut -d '_' -f1 | sed 's/^0//' | sed 's/^0//' | rev | cut -c 2- | rev`
pub_bperp=$LiCSAR_public/$tr/$frame/metadata/baselines

#regenerate baselines file if it exists:
#if [ -f baselines ]; then
# mv baselines baselines.bck
#fi

echo "preparing baselines file"
mk_bperp_file.sh

if [ ! -f $pub_bperp ]; then
 echo "generating baselines file in "$pub_bperp
 #echo "master    slave     bperp   btemp" > $pub_bperp
 cp baselines $pub_bperp
else
 echo "updating baselines file in "$pub_bperp
 if [ -f $pub_bperp.lock ]; then
  echo "error - seems the baselines are now in process of updating..cancelling for now"
  echo "you may check/delete file: "$pub_bperp".lock"
 else
  touch $pub_bperp.lock
  #tail -n+2 $pub_bperp > temp_bperp_file
  cp $pub_bperp temp_baselines
  #getting rid of potential nomaster error here:
  sed '/^0.0000/d' baselines >> temp_baselines
  #echo "master    slave     bperp   btemp" > $pub_bperp
  #sometimes there can be two/more different bperp lines. so 'sort -u' would not help, awk is very fast in here..
  awk -F" " '!_[$2]++' temp_baselines > $pub_bperp
  #sort -u temp_baselines > $pub_bperp
  rm temp_baselines
  rm $pub_bperp.lock
 fi
fi

#return the original bperp_file (in case someone wanted to use it)
#if [ -f baselines.bck ]; then
# mv baselines.bck baselines
#fi
