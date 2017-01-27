#!/bin/sh


rsyncToTmp() {
  [ -z "$1" ] && { echo "Usage: $0 index.gtf.bowtieIndex"; return; } 
  
  tmp=`mktemp -d`
   
  echo Start to rsync to tmp folder: $tmp
  for v in ${1//./ }; do
      echo working on $v
      #eval vx=\$$v
      vx=${!v}
      #echo real path is $vx

      if [ -e $vx ]; then
	#echo file or folder
	rsync -a ${vx%/} $tmp	  
      else
	#echo multiple files
	rsync -a $vx* $tmp
      fi      
      eval $v=$tmp/${vx##*/} 
    
   done

  echo rsync is done. Here is the tmp folder content: 
  myTmpDir=$tmp
  ls -l $myTmpDir
}
export -f rsyncToTmp

removeTmp() {
  #echo exit code is $?
  echo before delete tmp files
  ls -l $myTmpDir
  echo tmpdir is $myTmpDir
  rm -r $myTmpDir 
  echo after delete tmp files
  ls -l $myTmpDir  2>&1
  echo done
 
}	
export -f removeTmp
