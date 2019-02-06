#!/bin/bash
                
to=`cat ~/.forward`

flag=$1

# call by: sendJobFinishEmail fullPathFlag or sendJobFinishEmail out.txt err.txt script.txt 
[ -z "$2" ] && { out=$flag.out; script=$flag.sh; ss=$flag.success; } || { out="$1"; err="$2"; script="$3"; }

minimumsize=9000

actualsize=`wc -c $out`

[ -z "$2" ] && [ ! -f $ss ] && s="Subject: Failed: job id:$SLURM_JOBID name:$SLURM_JOB_NAME\n" ||  s="Subject: Success: job id:$SLURM_JOBID name:$SLURM_JOB_NAME\n"

[ ! -z "$2" ] && `grep jobSuccessfullyDone $out` &&   s="Subject: Success: job id:$SLURM_JOBID name:$SLURM_JOB_NAME\n" || s="Subject: Failed: job id:$SLURM_JOBID name:$SLURM_JOB_NAME\n" 

stat=`tail -n 1 $out`
[[ "$stat" == *COMPLETED* ]] && echo *Notice the sacct report above: while the main job is still running for sacct command, user task is completed. >> $out

if [ "${actualsize% *}" -ge "$minimumsize" ]; then
   toSend=`echo Job script content:; cat $script`
   toSend="$s\n$toSend\nOutput is too big for email. Please find output in: $out"  
   toSend="$toSend\n...\n`tail -n 6 $out`"
   [ -f "$err" ] && toSend="$toSend\n...Err:\n`tail -n 6 $err`"
else
   toSend=`echo Job script content:; cat $script; echo Job output:; cat $out`
   toSend="$s\n$toSend"
   [ -f "$err" ] && toSend="$toSend\n...Err:\n`tail -n 6 $err`"
fi

echo -e "$toSend" | sendmail $to 
