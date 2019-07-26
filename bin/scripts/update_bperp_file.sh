#/bin/bash
#to be run from frame directory!

frame=`pwd | rev | cut -d '/' -f1 | rev`
if [ ! `echo $frame | cut -d '_' -f3 | wc -c` -eq 7 ]; then echo "to be run from frame directory"; exit; fi

tr=`echo $frame | cut -d '_' -f1 | sed 's/^0//' | sed 's/^0//' | rev | cut -c 2- | rev`
pub_bperp=$LiCSAR_public/$tr/$frame/metadata/bperp_file

#regenerate bperp_file if it exists:
if [ -f bperp_file ]; then
 mv bperp_file bperp_file.bck
fi

echo "preparing bperp file"
mk_bperp_file.sh

if [ ! -f $pub_bperp ]; then
 echo "generating bperp_file in "$pub_bperp
 #echo "master    slave     bperp   btemp" > $pub_bperp
 cp bperp_file $pub_bperp
else
 echo "updating bperp_file in "$pub_bperp
 #tail -n+2 $pub_bperp > temp_bperp_file
 cp $pub_bperp temp_bperp_file
 cat bperp_file >> temp_bperp_file
 #echo "master    slave     bperp   btemp" > $pub_bperp
 sort -u temp_bperp_file >> $pub_bperp
 rm temp_bperp_file
fi

#return the original bperp_file (in case someone wanted to use it)
if [ -f bperp_file.bck ]; then
 mv bperp_file.bck bperp_file
fi
