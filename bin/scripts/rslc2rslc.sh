#! /bin/bash

# Find offset in azimuth between master and cropped master

master1=$1 #RSLC/${masterdate}/${masterdate} #$1
slave1=$2  #RSLC/${slavedate}/${masterdate} # $2
offset=$3 

# Size (Width/Length) of master
masterwid=`awk '$1 == "range_samples:" {print $2}' ${master1}.rslc.par`
masterlen=`awk '$1 == "azimuth_lines:" {print $2}' ${master1}.rslc.par`
# Size (Width/Length) of master
slavewid=`awk '$1 == "range_samples:" {print $2}' ${slave1}.rslc.par`
slavelen=`awk '$1 == "azimuth_lines:" {print $2}' ${slave1}.rslc.par`

# Python version
rslc2rslc.py ${master1}.rslc ${slave1}.rslc $offset 

# # Matlab version
# cat << _EOF_ > rslc2rslc.m
# slavecrop=freadbk('${slave1}.rslc',$slavelen,'cpxint16','b');
# slave=complex(zeros($masterlen,$masterwid));
# %slave(1+${offset}:size(slavecrop,1)+${offset},1+${offset}:size(slavecrop,2)+${offset})=slavecrop;
# slave(1+${offset}:size(slavecrop,1)+${offset},:)=slavecrop;
# fwritebk(slave,'${slave1}.full.rslc','cpxint16','b');
# _EOF_
# matlab -nojvm -nosplash < rslc2rslc.m

# Convert offset to time in the azimuth time of the slave!
# That's to adjust the start_time, center_time and end_time in the $slave.rslc.par file
az_line_time=`awk '$1 == "azimuth_line_time:" {print $2}' ${slave1}.rslc.par`  
slave_start_time=`awk '$1 == "start_time:" {print $2}' ${slave1}.rslc.par`  
slave_end_time=`awk '$1 == "end_time:" {print $2}' ${slave1}.rslc.par`  
start_time2=`echo $slave_start_time $offset $az_line_time | awk '{printf "%.6f", $1-($2*$3)}'` 
end_time2=`echo $start_time2 $masterlen $az_line_time | awk '{printf "%.6f", $1+($2*$3)}'` 
center_time2=`echo $end_time2 $start_time2 | awk '{printf "%.6f", $2+(($1-$2)/2)}'` 

awk '{if ($1=="range_samples:") print "range_samples:                 "'${masterwid}'""; else print $0}' ${slave1}.rslc.par > temp
awk '{if ($1=="azimuth_lines:") print "azimuth_lines:                 "'${masterlen}'""; else print $0}' temp > temp1
awk '{if ($1=="start_time:") printf "%s %.6f%s\n", "start_time:            ","'${start_time2}'","   s"; else print $0}' temp1 > temp2
awk '{if ($1=="center_time:") printf "%s %.6f%s\n", "center_time:           ","'${center_time2}'","   s"; else print $0}' temp2 > temp3
awk '{if ($1=="end_time:") printf "%s %.6f%s\n", "end_time:              ","'${end_time2}'","   s"; else print $0}' temp3 > temp4
#sed -i 's/lines_per_burst: *'$burst_smp'/lines_per_burst:   '$maxburst_smp'/g'  SLC/$yyyymmdd/$yyyymmddhhmmss*IW$n.slc.TOPS_par 
#sed -i 's/range_samples: *'$burst_smp'/lines_per_burst:   '$masterwid'/g'  ${slave1}.rslc.par 

# Rename products after azimuth shift
mv temp4 ${slave1}.rslc.par
mv temp1 ${slave1}.rslc.par
mv ${slave1}.full.rslc ${slave1}.rslc

rm -f temp*


