#!/bin/bash

##########################################################

#fix for bsub -> slurm
a=`which bsub 2>/dev/null`
if [ -z $a ]; then
 echo "setting bsub2slurm"
 alias bsub=bsub2slurm.sh
fi


function cdproc() {
    if [ -z $1 ]; then
        echo "please provide frame id";
        return 0
    fi
    frame=$1
    tr=`track_from_frame $frame`
    cd $LiCSAR_procdir/$tr/$frame
}
export -f cdproc


function cdpub() {
    if [ -z $1 ]; then
        echo "please provide frame id";
        return 0
    fi
    frame=$1
    tr=`track_from_frame $frame`
    cd $LiCSAR_public/$tr/$frame
}
export -f cdpub