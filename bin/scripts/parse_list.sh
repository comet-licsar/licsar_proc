#!/bin/bash

#################################################################
# spit out list of unique days contained within a list of S1 zip files
#
# Dependencies:
#
# Author: Jonathan Weiss
# Date: 2018/07/18
#################################################################

usage() { 
  echo " " 1>&2; 
  exit 1; 
}

if [ "$#" -ne 1 ]; then
  usage
else
  file=$1
fi

sed 's/\// /g' $file | awk '{print $10}' | cut -c 18-25 | sort | uniq
