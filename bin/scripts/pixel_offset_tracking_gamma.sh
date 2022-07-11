# how to mosaic stuff:

a=106D_05447_131313
b=106D_05646_131313

b=099A_05416_131313
a=099A_05615_131313
epoch1=20210520
epoch2=20210526

for epoch in $epoch1 $epoch2; do
    rm off.bagr 2>/dev/null
    mkdir RSLC/$epoch
    create_offset ../$a/RSLC/$epoch/$epoch.rslc.par ../$b/RSLC/$epoch/$epoch.rslc.par off.bagr 1 20 4 0
    SLC_cat ../$a/RSLC/$epoch/$epoch.rslc ../$b/RSLC/$epoch/$epoch.rslc ../$a/RSLC/$epoch/$epoch.rslc.par ../$b/RSLC/$epoch/$epoch.rslc.par off.bagr RSLC/$epoch/$epoch.rslc RSLC/$epoch/$epoch.rslc.par
    rm off.bagr 2>/dev/null
done





# merge them properly...:
path_frame=/home/home02/earmla/licsar_disk/chinaeq/desc
path_frame=/home/home02/earmla/licsar_disk/chinaeq/asc
cd $path_frame/mosaic
frame1=106D_05447_131313
frame2=106D_05646_131313

frame2=099A_05416_131313
frame1=099A_05615_131313

#try mosaicing through SLC_copy
# BAD OUTPUT!!!!!
for d in $epoch1 $epoch2; do
 for frame in $frame1 $frame2; do
  SLC_copy $path_frame/$frame/RSLC/$d/$d.rslc $path_frame/$frame/RSLC/$d/$d.rslc.par $path_frame/$frame/RSLC/$d/$d.rslc.scomplex $path_frame/$frame/RSLC/$d/$d.rslc.scomplex.par 3
 done
MLI_cat $path_frame/$frame1/RSLC/$d/$d.rslc.scomplex $path_frame/$frame2/RSLC/$d/$d.rslc.scomplex $path_frame/$frame1/RSLC/$d/$d.rslc.scomplex.par $path_frame/$frame2/RSLC/$d/$d.rslc.scomplex.par RSLC/$d/$d.rslc RSLC/$d/$d.rslc.par 1 0 0 0
done







# make pixel offset tracking

create_offset RSLC/$epoch1/$epoch1.rslc.par RSLC/$epoch2/$epoch2.rslc.par $epoch1'_'$epoch2.fullpix.off 1 1 1 0
#create_offset RSLC/$epoch1/$epoch1.rslc.par RSLC/$epoch2/$epoch2.rslc.par $epoch1'_'$epoch2.mlpix.off 1 20 4 0
time offset_pwr_tracking RSLC/$epoch1/$epoch1.rslc RSLC/$epoch2/$epoch2.rslc RSLC/$epoch1/$epoch1.rslc.par RSLC/$epoch2/$epoch2.rslc.par $epoch1'_'$epoch2.fullpix.off $epoch1'_'$epoch2.offsets $epoch1'_'$epoch2.corr
offset_tracking $epoch1'_'$epoch2.offsets $epoch1'_'$epoch2.corr RSLC/$epoch1/$epoch1.rslc.par $epoch1'_'$epoch2.fullpix.off $epoch1'_'$epoch2.disp_map  $epoch1'_'$epoch2.disp_val 1 - 0

# Extract Range offsets (stores as real and imaginary part)
widthoff=`grep range_samples $epoch1'_'$epoch2.fullpix.off | awk '{print $2}'`
cpx_to_real $epoch1'_'$epoch2.disp_map $epoch1'_'$epoch2.disp_map.rng $widthoff 0
swap_bytes $epoch1'_'$epoch2.disp_map.rng $epoch1'_'$epoch2.disp_map.rng.swapped 4

#now geocode it


#and make preview
cpxfiddle -w $widthoff -o sunraster -c jet -q normal -f r4 -M5/1 -r-2.5/2.5 $epoch1'_'$epoch2.disp_map.rng.swapped > $epoch1'_'$epoch2.disp_map.rng.ras
#display $epoch1'_'$epoch2.disp_map.rng.ras



cd $x
swap_bytes IFG/$epoch1'_'$epoch2/$epoch1'_'$epoch2.cc IFG/$epoch1'_'$epoch2/$epoch1'_'$epoch2.cc.swapped 4
swap_bytes IFG/$epoch1'_'$epoch2/$epoch1'_'$epoch2.filt.diff IFG/$epoch1'_'$epoch2/$epoch1'_'$epoch2.filt.diff.swapped 4
cd ..

3401
4622
for x in 106D_05447_131313  106D_05646_131313; do mkdir -p pokus/$x/IFG/$epoch1'_'$epoch2; cp $x/IFG/$epoch1'_'$epoch2/$epoch1'_'$epoch2.*swapped pokus/$x/IFG/$epoch1'_'$epoch2/.; cp $x/$epoch1'_'$epoch2.disp_map.rng.swapped pokus/$x/.; widthoff=`grep range_samples $x/$epoch1'_'$epoch2.fullpix.off | awk '{print $2}'`; echo $widthoff > pokus/$x/width.txt; done


#this is not needed anymore:

for d in $epoch1 $epoch2; do
 for frame in $frame1 $frame2; do
  rm $frame.$d.tab
  for iw in 1 2 3; do echo $path_frame/$frame/RSLC/$d/$d.IW$iw.rslc $path_frame/$frame/RSLC/$d/$d.IW$iw.rslc.par $path_frame/$frame/RSLC/$d/$d.IW$iw.rslc.TOPS_par >> $frame.$d.tab; done
 done
 rm $d.tab
 for iw in 1 2 3; do echo RSLC/$d/$d.IW$iw.rslc RSLC/$d/$d.IW$iw.rslc.par RSLC/$d/$d.IW$iw.rslc.TOPS_par >> $d.tab; done
 # merge them and generate mosaic
 SLC_cat_ScanSAR $frame1.$d.tab $frame2.$d.tab $d.tab
 SLC_mosaic_S1_TOPS $d.tab RSLC/$d/$d.rslc RSLC/$d/$d.rslc.par 20 4
done



