#!/bin/bash
#################################################################
# Configuration file to launch LiCSAR stack InSAR processing
#
# Dependencies: 
#   GAMMA 20151209 (http://www.gamma-rs.ch/)
#   GDAL (http://www.gdal.org/)
#   Python 2.7 (shipped via Canopy) - numpy and pyAPS
#
# Author: Pablo J. Gonzalez
# Date: 2016/03/14
#################################################################
# Load all the software tools for LiCSAR to work
source /nfs/a1/software/LiCSAR/LiCSAR_v1.0/LiCSAR_source.sh 
# Tell LiCSAR where the zip SAFE files are located and some general parameters
RAW_DIR=/nfs/a1/raw/sentinel/Taiwan/T069A
master=`cat master_date.txt` # Force a given date to be the master, otherwise leave blank, and it uses the first date available
rlks=5                      # Multilooking factor in range 
azlks=1                      # Multilooking factor in azimuth 
outres=0.0002                 # Spatial resolution in degrees of the generated products [e.g., 20:4 = 0.001 = 100 metres / 5:1 = 0.0002 = 20 metres]
specdiv=1                    # Use spectral diversity once to improve coregistration
winsize=64                   # Conservative window size
geo_proc=0                   # Process the interferograms for filtering and unwrapping in geographical coordinates
aps_correction=0             # If you want to run atmospheric correction flag=1 (only for scenes older than 2 months)!!!
polygon_file=Taiwan.xy       # File containing the burst edges or a rectangular region to be processed
burstid_file=$(basename $polygon_file .xy)_burst_ids.txt
starting_step=1              # Execute controlled steps: 1 to 6 means from zip files to unfiltered interferograms
ending_step=10               # Execute controlled steps: 1 to 6 means from zip files to unfiltered interferograms

# Number of parallel threads the user wants to run
export OMP_NUM_THREADS=8         

# The first step is to check that either you introduced a valid master or polygon_file
#echo "init_LiCSAR -m $master -p $polygon_file "
plot_bursts.gmt ${RAW_DIR} $master $polygon_file $burstid_file
ls ${RAW_DIR}/S1A_*.zip > zipfile.list

#echo "LiCSAR -i zipfile.list -s $starting_step -e $ending_step -m $master -r $rlks -a $azlks -o $outres -d $specdiv -w $winsize -z $aps_correction -g $geo_proc -f $polygon_file "
LiCSAR -i zipfile.list -s $starting_step -e $ending_step -m $master -r $rlks -a $azlks -o $outres -d $specdiv -w $winsize -z $aps_correction -g $geo_proc -f $polygon_file

