## The script is used to generate the required input file for the GACOS.sh API
## It takes the frame id as the input
# optionaly, you may add second parameter to input and that will be a text file with dates to process
# e.g. 
# 20190401
# 20220202 ...
frame=$1
input=''
if [ ! -z $2 ]; then
 if [ ! -f $2 ]; then echo "if you use 2nd parameter, it should be an existing txt file with a list of epochs"; exit; fi
 input=`realpath $2`
fi

source $LiCSARpath/lib/LiCSAR_bash_lib.sh
# request gacos only for epochs older than $minbtemp
# update 2023: GACOS data are available already in 3 days lag only
minbtemp=3

prod_dir=$LiCSAR_public
work_dir="/gws/nopw/j04/nceo_geohazards_vol2/LiCS/temp/GACOS"

if [ -z $prod_dir ]; then
 echo "you probably did not load licsar_proc or licsar_framebatch"
 exit
fi

cd $work_dir
rm $frame.inp 2>/dev/null

touch $frame.inp
echo "$work_dir/$frame" >> $work_dir/$frame.inp

track=`echo $frame | cut -c -3 | sed 's/^0//' | sed 's/^0//'`


if [ -f $prod_dir/$track/$frame/metadata/metadata.txt ]; then
   source $prod_dir/$track/$frame/metadata/metadata.txt
   master=$(echo $master) 
else
   echo "File $FILE does not exist."
   echo $prod_dir/$track/$frame/metadata/metadata.txt
   rm $frame.inp
   exit 1
fi

hgttif=$prod_dir/$track/$frame/metadata/$frame'.geo.hgt.tif'
if [ -f $hgttif ]; then
   gmt grdinfo -C -M $hgttif | gawk {'print $4" "$5" "$2" "$3'} >> $work_dir/$frame.inp
   #minlat=`awk 'BEGIN{a=1000}{if ($2<0+a) a=$2} END{print a}' $prod_dir/$track/$frame/metadata/$frame'-poly.txt'`
   #maxlat=`awk 'BEGIN{a=   0}{if ($2>0+a) a=$2} END{print a}' $prod_dir/$track/$frame/metadata/$frame'-poly.txt'`
   #minlong=`awk 'BEGIN{a=1000}{if ($1<0+a) a=$1} END{print a}' $prod_dir/$track/$frame/metadata/$frame'-poly.txt'`
   #maxlong=`awk 'BEGIN{a=   0}{if ($1>0+a) a=$1} END{print a}' $prod_dir/$track/$frame/metadata/$frame'-poly.txt'`
   #echo $minlat $maxlat $minlong $maxlong >> $work_dir/$frame.inp
else
   echo "File $FILE does not exist."
   echo $hgttif
   rm $frame.inp
   exit 1
fi

if [ -f $prod_dir/$track/$frame/metadata/metadata.txt ]; then
   source $prod_dir/$track/$frame/metadata/metadata.txt
   h=$(echo $center_time | cut -d: -f1)
   m=$(echo $center_time | cut -d: -f2)
   s=$(echo $center_time | cut -d: -f3)
   echo $h+ $m/60+ $s/3600 | bc -l >>$frame.inp 
else
   echo "File $FILE does not exist."
   echo $prod_dir/$track/$frame/metadata/metadata.txt
   rm $frame.inp
   exit 1
fi

if [ ! -z $input ]; then
 cp $input $frame.inp.tmp2;
else
 ls $prod_dir/$track/$frame/int*/20*_* -d | rev | cut -d '/' -f1 | rev | cut -d '_' -f1 > $frame.inp.tmp
 ls $prod_dir/$track/$frame/int*/20*_* -d | rev | cut -d '/' -f1 | rev | cut -d '_' -f2 >> $frame.inp.tmp
 ls $prod_dir/$track/$frame/epochs/20* -d | rev | cut -d '/' -f1 | rev >> $frame.inp.tmp
 sort -u $frame.inp.tmp >$frame.inp.tmp2
fi
for epoch in `cat $frame.inp.tmp2`; do
    if [ ! -f $prod_dir/$track/$frame/epochs/$epoch/$epoch.sltd.geo.tif ]; then
     if [ `datediff $epoch` -ge $minbtemp ]; then
      echo $epoch >> $frame.inp
     fi
    fi
done
rm $frame.inp.tmp $frame.inp.tmp2 2>/dev/null

exit



if [ -f $prod_dir/$track/$frame/metadata/baselines ]; then
   for epoch in `awk '{print $2}' $prod_dir/$track/$frame/metadata/baselines`; do
    if [ ! -f $prod_dir/$track/$frame/epochs/$epoch/$epoch.sltd.geo.tif ]; then
     echo $epoch >> $frame.inp
    fi
   done
else
   echo "File $FILE does not exist."
   echo $prod_dir/$track/$frame/metadata/baselines
   rm $frame.inp
   exit 1
fi


