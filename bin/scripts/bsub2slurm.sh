#!/bin/bash
# saner programming env: these switches turn some bugs into errors
set -o errexit -o pipefail -o noclobber -o nounset


# IMPORTANT TODO:
# use jid1=$(sbatch jobscript | tr -dc '0-9')  to get job id for waiting script


# -allow a command to fail with !’s side effect on errexit
# -use return value from ${PIPESTATUS[0]}, because ! hosed $?
! getopt --test > /dev/null 
#if [[ ${PIPESTATUS[0]} -ne 4 ]]; then
#    echo 'I’m sorry, `getopt --test` failed in this environment.'
#    exit 1
#fi

#OPTIONS=dfo:v
#LONGOPTS=debug,force,output:,verbose

# -regarding ! and PIPESTATUS see above
# -temporarily store output to be able to check for errors
# -activate quoting/enhanced mode (e.g. by writing out “--options”)
# -pass arguments only via   -- "$@"   to separate them correctly
#! PARSED=$(getopt --options=$OPTIONS --longoptions=$LONGOPTS --name "$0" -- "$@")
#if [[ ${PIPESTATUS[0]} -ne 0 ]]; then
    # e.g. return value is 1
    #  then getopt has complained about wrong arguments to stdout
#    exit 2
#fi
# read getopt’s output this way to handle the quoting right:
#eval set -- "$PARSED"

#d=n f=n v=n outFile=-
cmd='sbatch '
# now enjoy the options in order and nicely split until we see --
addedextrapost=''
addedextrapre=''
while true; do
    case "$1" in
        -q)
            cmd=$cmd' -p '"$2"
            if [ $2 == 'cpom-comet' ]; then
             cmd=$cmd' --account=cpom-comet'
            fi
            shift 2
            ;;
        -W)
            cmd=$cmd' --time='"$2"':00'
            shift 2
            ;;
        -J)
            cmd=$cmd' --job-name='"$2"
            shift 2
            ;;
        -o)
            cmd=$cmd' --output='"$2"
            shift 2
            ;;
        -e)
            cmd=$cmd' --error='"$2"
            shift 2
            ;;
        -n)
            cmd=$cmd' -n '"$2"
            shift 2
            ;;
        -E)
            cmd=$cmd' --prolog '\"$2\"
            addedextrapre=$2'; '
            shift 2
            ;;
        -Ep)
            #cmd=$cmd' --epilog '\"$2\"
            addedextrapost='; '$2
            shift 2
            ;;
        -R)
            mem=`echo $2 | cut -d '=' -f2 | cut -d ']' -f1`
            cmd=$cmd' --mem='$mem
            shift 2
            ;;
        -M)
            cmd=$cmd' --mem='$2
            shift 2
            ;;
        -w)
            #e.g. -w "ended(framebatch_02_coreg_559619) && ended(framebatch_02_coreg_559621)"
            vars=$2
            jobids=''
            for myJOBNAME in `echo $vars | sed 's/ended(//g' | sed 's/)//g' | sed 's/\&\&//g' | tr "'" " "`; do
             #this way the jobid can be really 'historic'
             #jobid=$(sacct -n --format="JobID" --name $myJOBNAME | head -n1 | cut -d '.' -f1)
             jobid=$(squeue -n $myJOBNAME | sed '/JOBID/d' | head -n1 | gawk {'print $1'})
             #jobid=$myJOBNAME
             if [ ! -z $jobid ]; then
               jobids=$jobids':'$jobid
             else
               #echo "ERROR: dependency not satisfied - seems job "$myJOBNAME" is not active.."
               #echo "trying with archived processing info, but expect problems"
               jobid=$(sacct -n --name $myJOBNAME | head -n1 | gawk {'print $1'})
               if [ ! -z $jobid ]; then
                jobids=$jobids':'$jobid
               fi
               #exit
             fi
            done
            cmd=$cmd' --dependency=afterany'$jobids
            shift 2
            ;;
        -v|--verbose)
            v=y
            shift
            ;;
        *)
            cmd=`echo $cmd' --wrap='\"$addedextrapre $* $addedextrapost\"`
            #echo 'tu'
            #exit 3
            shift
            break
            ;;
    esac
done

#echo $cmd
eval $cmd

# handle non-option arguments
##if [[ $# -ne 1 ]]; then
#    echo "$0: A single input file is required."
#    exit 4
#fi


