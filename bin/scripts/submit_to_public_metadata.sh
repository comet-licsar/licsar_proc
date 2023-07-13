#!/bin/bash

if [ -z $1 ]; then echo "parameter is FRAME"; exit; fi

frame=$1
track=`echo $frame | cut -c -3 | sed 's/^0//' | sed 's/^0//'`
outfile=$LiCSAR_public/$track/$frame/metadata/metadata.txt

master=`ls $LiCSAR_procdir/$track/$frame/geo/*[0-9].hgt | rev | cut -d '/' -f1 | rev | cut -d '.' -f1 | head -n1`

#write master time to the outfile
echo "master="$master > $outfile

#write coarse acquisition time (e.g. for GACOS purposes)
time=`grep center_time $LiCSAR_procdir/$track/$frame/SLC/$master/$master.slc.par | gawk {'print $2'}`
time_hours=`python -c "import datetime as dt; time=dt.timedelta(seconds="$time"); print(str(time))"`
echo "center_time="$time_hours >> $outfile

#write heading - e.g. for evaluation of SD values
#hfile=`ls $track/$frame/geo/*.heading.grd 2>/dev/null`
hfile=`ls $LiCSAR_procdir/$track/$frame/geo/*.geo.deg.heading 2>/dev/null`
if [ -z $hfile ]; then 
 echo "frame "$frame" has no heading file";
else
 heading=`python3 -c "import numpy as np; a = np.fromfile('"$hfile"', dtype=np.float32); print(a.byteswap().mean())"`
 #heading=`gmt grdinfo -L2 $hfile | grep mean | gawk {'print $3'}`
 echo "heading="$heading >> $outfile
fi

#write incidence angle
inc_angle=`grep incidence_angle $LiCSAR_procdir/$track/$frame/SLC/$master/$master.slc.par | gawk {'print $2'}`
echo "avg_incidence_angle="$inc_angle >> $outfile

#write resolution in azimuth, and range
valuee=`grep azimuth_pixel_spacing $LiCSAR_procdir/$track/$frame/SLC/$master/$master.slc.par | gawk {'print $2'}`
echo "azimuth_resolution="$valuee >> $outfile
valuee=`grep range_pixel_spacing $LiCSAR_procdir/$track/$frame/SLC/$master/$master.slc.par | gawk {'print $2'}`
echo "range_resolution="$valuee >> $outfile

#write center_range_slc
valuee=`grep center_range_slc $LiCSAR_procdir/$track/$frame/SLC/$master/$master.slc.par | gawk {'print $2'}`
echo "centre_range_m="$valuee >> $outfile

#write average height:
heifile=`ls $LiCSAR_procdir/$track/$frame/geo/20??????.hgt 2>/dev/null`
if [ -z $heifile ]; then 
 echo "frame "$frame" has no height file - error";
else
 hei=`python3 -c "import numpy as np; a = np.fromfile('"$heifile"', dtype=np.float32); print(a.byteswap().mean())"`
 echo "avg_height="$hei >> $outfile
fi

# write info on DEM used here
# cd $LiCSAR_public; for track in `ls`; do for frame in `ls $track`; do
# outfile=$track/$frame/metadata/metadata.txt; if [ -f $outfile ]; then
demlog=`ls $LiCSAR_procdir/$track/$frame/01_doDEMcrop.log 2>/dev/null`
if [ -z $demlog ]; then
  echo "frame "$frame" has no DEM log file - error";
else
  if [ `grep -c Copernicus $demlog` -gt 0 ]; then demstr="CopernicusDEM_30m";
  elif [ `grep -c TanDEM $demlog` -gt 0 ]; then demstr="TanDEMX_90m";
  elif [ `grep -c tdm2insar $demlog` -gt 0 ]; then demstr="TanDEMX_90m";
  elif [ `grep -c gdem2insar $demlog` -gt 0 ]; then demstr="AsterGDEM_30m";
  elif [ `grep -c srtm2insar $demlog` -gt 0 ]; then demstr="SRTM_30m";
  fi
  echo "applied_DEM="$demstr >> $outfile
fi
# fi; done; done
