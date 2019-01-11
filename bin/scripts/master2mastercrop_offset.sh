#! /bin/bash
###########################################################
# Find offset in azimuth between master and cropped master
#
# Knowing the azimuth line time rate and time difference,
# We can guess the number of lines to displace the cropped
# images into the full master image
#
# This assumes that the width of the files are the same
#
# Author: Pablo J. Gonzalez
# Date: 2016/04/06
###########################################################

master=$1 #RSLC/${masterdate}/${masterdate} #$1
extmaster=$2 # rslc
master_cropped=$3 #RSLC/${slavedate}/${masterdate} # $2
extmaster_cropped=$4 # slc or rslc

# determine lines offset between start of burst1 and start of burst2
az_start_time1=`awk '$1 == "start_time:" {print $2}' ${master}.${extmaster}.par`      
az_start_time2=`awk '$1 == "start_time:" {print $2}' ${master_cropped}.${extmaster_cropped}.par`        
az_line_time=`awk '$1 == "azimuth_line_time:" {print $2}' ${master}.${extmaster}.par`      

#echo "$az_start_time1 $az_start_time2 $az_line_time"
#lines_offset_float=`echo "$az_start_time1 $ax_start_time2 $azimuth_line_time" | awk '{printf "%f", (($2-$1)/$3)}'`
lines_offset=`echo "$az_start_time1 $az_start_time2 $az_line_time" | awk '{printf "%d", (0.5+($2-$1)/$3)}'`

echo $lines_offset

  