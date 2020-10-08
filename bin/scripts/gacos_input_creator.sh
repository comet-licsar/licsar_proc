## The script is used to generate the required input file for the GACOS.sh API
## It takes the frame id as the input

frame=$1

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

if [ -f $prod_dir/$track/$frame/metadata/baselines ]; then
   for epoch in `awk '{print $2}' $prod_dir/$track/$frame/metadata/baselines`; do
    if [ ! -f $prod_dir/$track/$frame/epochs/$epoch/$epoch.ztd.geo.tif ]; then
     echo $epoch >> $frame.inp
    fi
   done
else
   echo "File $FILE does not exist."
   echo $prod_dir/$track/$frame/metadata/baselines
   rm $frame.inp
   exit 1
fi


