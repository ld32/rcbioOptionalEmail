#!/bin/sh

Usage="Usage: $0 full_path_to_flag_folder  flag(job name or job id) [remove_success_flag_file]  \n  Note: this script will go through job id list file, find the downstream jobs, and kill them if the current job failed or remove success flag file. "

path=$1

[ -z "$2" ] && { echo -e "job name or job id  can not be empty! \n$Usage"; exit; }
        
           # kill downsteam jobs                          or  remove .successflag 
[ -z "$3" ] && { kill=""; text=`cat $1/alljobs.jid`; }  || { kill=kill; text=`cat $1/alljobs.jid.first`; }

[ -f $path/alljobs.jid ] || { echo -e "job id file $path/alljobs.jid does not exist\n$Usage" >> $path/$2.out; exit; }

[ -f $path/$2.success ] && echo kill jobs should not run, because job has finished successfully. >> $path/$2.out && exit 
      
job=$2

echo checking to see if job finished successfully. job name is $job >> $path/$job.out

IFS=$' ';  

# check the third column for the job name, then find the the job id in column 1
id=`echo $text | awk '{if ($3 ~ /'"$job/"') print $1,$3;}'`

# job name can not find, it maybe a job id  instead of job name 
[ -z "$id" ] && id=`echo $text | awk '{if ($1 ~ /'"$job/"') print $1,$3;}'`
#echo job id $id

[ -z "$id" ] && { echo -e "job name or job id  can not be empty! \n$Usage"; exit; }

declare -A seen

function findID {
    #echo function start for $1
    ids="$1" 
    [ -z "$ids" ] &&  return

    IFS=$'\n';   
    for k in $ids; do
    	i=${k% *}; j=${k#* };
        #echo working on $k ,   get  i = $i , j = $j 
        if [ ! -z  $itself ]; then 
           [ -z "$kill" ] && echo killing downstream job or remove success flag file : $i $j 

            [ -z "$kill" ] && echo killing downstream job or remove success flag file : $i $j  >> $path/$job.out 
            if [ ! "${seen[$i]}" ]; then
              seen[$i]=1
              [ -z "$kill" ] && bkill $i >> $path/$job.out 2>&1 && touch $path/$j.killed || { echo rm $path/$j.success; rm $path/$j.success 2>/dev/null; } 
            
            else
              echo job killed before >> $path/$job.out
              echo job killed before
            fi
        fi
        itself=false  
       
        IFS=$' ';  
    	ids=`echo $text | awk '{if ($2 ~ /'"$i/"') print $1,$3;}'`
        #echo working on $i find: $ids | tr '\n' ' '; echo

        findID "$ids"
    done
}

findID "$id"

