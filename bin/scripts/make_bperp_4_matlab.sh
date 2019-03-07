#!/bin/bash

# create perp baseline file for plotting in matlab using v long temporal/perp baselines
# searches for the successfully unwrapped ifgs and cross-checks with above
# outputs new perp baseline file (location_frame) for plotting in matlab
#
# using plot_sb_list_gamma.m from Pablo to plot the output file
#

dir=`pwd`
frame=`basename ${dir}`
#location=`echo ${dir} | sed 's/\// /g' | awk '{print $(NF-1)}'`

#################################################################
# Dependencies:
# base_calc (GAMMA routine)

# Author: Jonathan Weiss
# Date: July 2018
#################################################################

echo working on frame $frame

#this bit actually calls the gamma base_calc routine
rm bperp* $frame'_bp.list' ifgs_done* SLC_* base_calc* 00_doB* 2>/dev/null

mk_bperp_jw.sh 5000 600
awk '{print $1,$2"_"$3,$4,$5,$6,$7,$8,$9}' bperp_file > bperp_file1

rm ifgs_done.list 2>/dev/null
#ls IFG/*/*unw.ras | sed 's/\// /g' | awk '{print $2}' > ifgs_done.list
ls IFG/*/*.unw | sed 's/\// /g' | awk '{print $2}' > ifgs_done.list

touch ${frame}_bp.list
touch temp
for ifg in `cat ifgs_done.list`; do grep $ifg bperp_file1 >> temp; done
sed 's/_/ /g' temp >> ${frame}_bp.list
\rm temp

done=`wc -l ifgs_done.list | awk '{print $1}'`
max=`wc -l ${frame}_bp.list | awk '{print $1}'`
tot=`wc -l bperp_file1 | awk '{print $1}'`


echo ${done} ifgs out of ${max} ifgs were found in bperp_file "(n_total=$tot)"

