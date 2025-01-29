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
v=0
lotusversion=1
if [ `hostname` == 'host839.jc.rl.ac.uk' ]; then
  echo "using LOTUS2"
  lotusversion=2
  qos='standard'
  hours=0
fi

while true; do
    case "$1" in
        -q)
            if [ $lotusversion -gt 1 ]; then
              if [ $2 == 'cpom-comet' ] || [ $2 == 'comet' ] || [ $2 == 'comet_lics' ] || [ $2 == 'comet_responder' ]; then
                cmd=$cmd' --account=comet_lics'
              else
                cmd=$cmd' --account=nceo_geohazards'
              fi
            else
              cmd=$cmd' -p '"$2"
              if [ $2 == 'cpom-comet' ] || [ $2 == 'comet' ] || [ $2 == 'comet_lics' ]; then
  #             cmd=$cmd' --account=cpom-comet'
                cmd=$cmd' --account=comet --partition=comet'
              fi
              if [ $2 == 'comet_responder' ]; then
               cmd=$cmd' --account=comet --partition=comet_responder'
              fi
            fi
            shift 2
            ;;
        -W)
            cmd=$cmd' --time='"$2"':00'
            hours=$2
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
            #cmd=$cmd' -n '"$2"
            cmd=$cmd' --cpus-per-task='"$2"
            if [ $2 -gt 1 ]; then
              qos='highres'
            fi
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
             jobid=''
             count=1
             max_count=1  # 2024 - bottleneck. doing it only once
             jobid=$(squeue -h --name=${myJOBNAME} --format='%i' | tail -n1)
             while [ -z "${jobid}" ] && [ ${count} -le ${max_count} ]
             do
                #echo "job not found, trying again - attempt "$count"/"$max_count
                #echo "ID of job is: "$myJOBNAME
                echo "trying by: squeue -h --format='%i' --name="$myJOBNAME
                # great solution by Rich Rigby! as sometimes jobs were not found...
                jobid=$(squeue -h --name=${myJOBNAME} --format='%i' | tail -n1)
                count=$((${count}+1))
                #sleep 3
             done
             #jobid=$(squeue -n $myJOBNAME | sed '/JOBID/d' | head -n1 | gawk {'print $1'})
             #jobid=$myJOBNAME
             if [ ! -z $jobid ]; then
               jobids=$jobids':'$jobid
             else
               echo "Job ID not found"
               #echo "ERROR: dependency not satisfied - seems job "$myJOBNAME" is not active.."
               #echo "trying with archived processing info, but expect problems"
               #jobid=$(sacct -n --name $myJOBNAME | head -n1 | gawk {'print $1'})
               domore=''
               if [ ! -z $domore ]; then
                 jobid=$(sacct -n --format=jobid --name=${myJOBNAME} | egrep '^[0-9]+\s' | sort -n | tail -n 1 | sed 's/ //g')
                 echo "trying alternative solution to figure jobID using:"
                 echo "sacct -n --format=jobid --name="$myJOBNAME" | egrep '^[0-9]+\s' | sort -n | tail -n 1 | sed 's/ //g'"
                 if [ ! -z $jobid ]; then
                  jobids=$jobids':'$jobid
                 else
                  echo "Job still not found - skipping this dependency"
                 fi
               fi
               #exit
             fi
            done
            if [ ! -z $jobids ]; then
             cmd=$cmd' --kill-on-invalid-dep=no --dependency=afterany'$jobids
            else
             echo "no jobs identified for dependency, sending to start without waiting"
            fi
            shift 2
            ;;
        -v|--verbose)
            v=1
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

if [ $lotusversion -gt 1 ]; then
  if [ $qos == 'standard' ]; then
    if [ $hours -gt 47 ]; then
      qos='long';
    elif [ $hours -gt 23 ]; then
      qos='highres'
    elif [ $hours -lt 5 ]; then
      qos='short'
    fi
  fi
  cmd=$cmd' --partition='$qos' --qos='$qos
fi

if [ $v == 1 ]; then
 echo $cmd
fi

eval $cmd

# handle non-option arguments
##if [[ $# -ne 1 ]]; then
#    echo "$0: A single input file is required."
#    exit 4
#fi


