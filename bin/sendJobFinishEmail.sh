#!/bin/bash
                
to=`cat ~/.forward`

flag=$1

minimumsize=9000

actualsize=`wc -c $flag.out`

[ -f $flag.failed ] && s="Subject: Failed: job id:$SLURM_JOBID name:$SLURM_JOB_NAME\n" ||  s="Subject: Success: job id:$SLURM_JOBID name:$SLURM_JOB_NAME\n"

stat=`tail -n 1 $flag.out`
[[ "$stat" == *COMPLETED* ]] && echo *Notice the sacct report above: while the main job is still running for sacct command, user task is completed. >> $flag.out

if [ "${actualsize% *}" -ge "$minimumsize" ]; then
   toSend=`echo Job script content:; cat $flag.sh`
   toSend="$s\n$toSend\nOutput is too big for email. Please find output in: $flag.out"  
else
   toSend=`echo Job script content:; cat $flag.sh; echo Job output:; cat $flag.out`
   toSend="$s\n$toSend"
fi

echo -e "$toSend" | sendmail $to 
