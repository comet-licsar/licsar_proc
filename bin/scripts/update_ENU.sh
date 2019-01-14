PATH1=/gws/nopw/j04/nceo_geohazards_vol1/projects/LiCS/proc/current
PATH_PUBLIC=/gws/nopw/j04/nceo_geohazards_vol1/public/LiCSAR_products

#update all
for TRACK in `ls $PATH_PUBLIC/* -d | rev | cut -d '/' -f1 | rev`; do
 for FRAME in `ls $PATH_PUBLIC/$TRACK/* -d | rev | cut -d '/' -f1 | rev`; do
  submit_lookangles.py -f $FRAME -t $TRACK > $FRAME'.log' 2> $FRAME'.err';
  echo $FRAME;
 done;
done

#correction of problems:
for FRAME in `grep 'cannot open file' *err | cut -d ':' -f1 | uniq | cut -d '.' -f1`; do
 echo "There was error in frame "$FRAME". Reprocessing"
 TRACK=`echo $FRAME | cut -d '_' -f1 | rev | cut -c 2- | rev | sed 's/^0//' | sed 's/^0//'`
 cd $PATH1/$TRACK/$FRAME
 MASTER=`ls geo/20*.lt | cut -d '/' -f2 | cut -d '.' -f1`
 cd $PATH1/$TRACK/$FRAME/RSLC/$MASTER
 for x in `ls ../../SLC/$MASTER/*`; do ln -s $x; done
 for x in `ls *slc*`; do mv $x `echo $x | sed 's/slc/rslc/'`; done
 submit_lookangles.py -f $FRAME -t $TRACK > $FRAME'.log' 2> $FRAME'.err'
done

#changing names
for track in `ls $PATH_PUBLIC/* -d | rev | cut -d '/' -f1 | rev`; do
 cd $PATH_PUBLIC/$track;
 for frame in `ls *_*_* -d`; do
  for tifpath in `ls $frame/metadata/*.geo.*.tif 2>/dev/null`; do
   nametif=`echo $tifpath | rev | cut -d '/' -f1 | rev`;
   master=`echo $nametif | cut -d '.' -f1`;
   if [ ! $master == $frame ]; then
    mv -f $tifpath $frame/metadata/`echo $nametif | sed 's/'$master'/'$frame'/'`;
   fi
  done;
 done;
 cd ..;
done

