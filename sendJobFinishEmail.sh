#!/bin/bash
                
to=`cat ~/.forward`

flag=$1

minimumsize=9000

actualsize=`wc -c $flag.out`

[ -f $flag.failed ] && s="Subject: Failed: id:$SLURM_JOBID name:$SLURM_JOB_NAME\n" ||  s="Subject: Success: id:$SLURM_JOBID name:$SLURM_JOB_NAME\n"

if [ "${actualsize% *}" -ge "$minimumsize" ]; then
   toSend=`cat $flag.sh`
   toSend="$s\n$toSend\nOutput is too big for email. Please find output in: $flag.out"  
else
   toSend=`cat $flag.sh $flag.out`
   toSend="$s\n$toSend"
fi

echo -e "$toSend" | sendmail $to 