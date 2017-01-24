#!/bin/sh

rsyncToTmp() {
  [ -z "$1" ] && { echo "Usage: $0 index.gtf.bowtieIndex"; return; } 
  
  for v in ${1//./ }; do
      echo working on $v
      #eval vx=\$$v
      vx=`realpath ${!v}` 
      myTmpDir=/tmp/${vx%/*}
      mkdir -p $myTmpDir
      if [[ $vx == */ ]]; then
         noEndingSlash=${vx%/}
         rsync -a $noEndingSlash  $myTmpDir
         eval $v=$myTmpDir/${noEndingSlash##*/}/
      else
         rsync -a ${vx}*  $myTmpDir
         eval $v=$myTmpDir/${vx##*/}         
      fi
      echo rsync done, from $vx to $myTmpDir
      echo Here is the tmp folder content:
      ls -l $myTmpDir
  done

}
