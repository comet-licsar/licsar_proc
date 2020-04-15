#!/bin/bash

if [ $# -lt 1 ]; then
  echo "----------------------------------------------------------------------"
  echo "GAOCS API Script, to run gacos on your local machine"
  echo "Contact chen.yu@ncl.ac.uk for any inquiries"
  echo ""
  echo "Note : Change username and password to your own"
  echo "Note : Must install curl"
  echo ""
  echo "Usage:"
  echo "      gacosapi.bash inputfile "
  echo "      inputfile is formated as a text file"
  echo "      Line 1: result file name (string)"
  echo "      Line 2: minlat maxlat minlon maxlon (float)"
  echo "      Line 3: time of the day in hours (float)"
  echo "      Line 4-end: date list (one date per line, formated as YYYYMMDD) "
  echo "An input file example:"
  echo "      mygacoscorrection"
  echo "      10.1 10.9 101.2 102.3"
  echo "      10.5"
  echo "      20180105"
  echo "      20180117"
  echo "      20180129"
  echo "-----------------------------------------------------------------------"
  exit
fi

username=gacos01
keyfile=$LiCSAR_configpath/id_rsa_gacos
inputfile=$1
pro_time=`date +%Y%m%d%H%M%S`


chmod 500 $keyfile
resultfile=`sed -n '1,1p' $inputfile `
col=`echo $resultfile | awk '{print NF}'`
if [ $col -ne 1 ]; then
   echo "ERROR - INFO - $timestr - RESULT FILE NAME WRONG"
   exit
fi

systemfile=API-${username}-${pro_time}.txt

outname_base=`sed -n '1,1p' $inputfile `
outname=`basename ${outname_base}`
outdir=`dirname ${outname_base}`
if [ ! -e $outdir ]; then 
   mkdir -p $outdir
fi
resultfile=${outname}_$pro_time
info=`sed -n '2,1p' $inputfile `
minlat=`echo $info | awk '{print $1}'`
maxlat=`echo $info | awk '{print $2}'`
minlon=`echo $info | awk '{print $3}'`
maxlon=`echo $info | awk '{print $4}'`
timeofday=`sed -n '3,1p' $inputfile `
echo $resultfile > $systemfile
echo $minlat $maxlat $minlon $maxlon >> $systemfile
echo $timeofday >> $systemfile



timestr=`date +%Y-%m-%d-%H:%M:%S`
echo "GACOSAPI - INFO - $timestr - GACOS Client version: 1.5"
echo "GACOSAPI - INFO - $timestr - Welcome user ${username}"
echo "GACOSAPI - INFO - $timestr - Input filename   = ${inputfile}"
echo "GACOSAPI - INFO - $timestr - Result filename  = ${outname_base}.tar.gz"
echo "GACOSAPI - INFO - $timestr - AREA             = ${minlat}/${maxlat}/${minlon}/${maxlon}"
echo "GACOSAPI - INFO - $timestr - TIME             = ${timeofday}"
echo "GACOSAPI - INFO - $timestr - DATELIST         :"
i=0
cat $inputfile | while read line
do
   i=`echo $i | awk '{print$1+1}'`
   if [ $i -le 3 ]; then continue; fi
   echo "        ${line}"
   echo $line >> $systemfile
done
timestr=`date +%Y-%m-%d-%H:%M:%S`
echo "GACOSAPI - INFO - $timestr - Request submitted"

echo finished >> ${systemfile}.finish


sftp -oIdentityFile=${keyfile} ${username}@sage-gacos-transfer.ncl.ac.uk:/${username}/inward<< ! > /dev/null 2>&1
put $systemfile
put ${systemfile}.finish
chmod 755 $systemfile
chmod 755 ${systemfile}.finish
bye
!

rm -rf ${systemfile}.finish ${systemfile}

resultfile=${resultfile}.tar.gz

sys_time=`date +%Y-%m-%d-%H:%M:%S`
echo "GACOSAPI - INFO - Start processing at $sys_time         "
while true
do
if [ -e ${resultfile}.finish ]; then

sftp -oIdentityFile=${keyfile} ${username}@sage-gacos-transfer.ncl.ac.uk:/${username}/outward<<EOT
get $resultfile
bye
EOT

break
fi
sleep 10s
sftp -oIdentityFile=${keyfile} ${username}@sage-gacos-transfer.ncl.ac.uk:/${username}/outward<< ! > /dev/null 2>&1
get ${resultfile}.finish
bye
!

done

rm ${resultfile}.finish
mv $resultfile ${outname_base}.tar.gz
sys_time=`date +%Y-%m-%d-%H:%M:%S`
echo "GACOSAPI - INFO - Finish processing at $sys_time         "
