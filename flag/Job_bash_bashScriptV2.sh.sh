#!/bin/bash
#Commands:
touch /home/ld32/rcbio/flag/Job_bash_bashScriptV2.sh.start
srun -n 1 bash -e -c "{ set -e; bash bashScriptV2.sh; } && touch /home/ld32/rcbio/flag/Job_bash_bashScriptV2.sh.success || touch /home/ld32/rcbio/flag/Job_bash_bashScriptV2.sh.failed"
sleep 5
echo Job done. Summary:
sacct --format=JobID,Submit,Start,End,State,Partition,ReqTRES%20,Timelimit,CPUTime,MaxRSS,NodeList%30 --units=M -j $SLURM_JOBID
sendJobFinishEmail.sh  /home/ld32/rcbio/flag/Job_bash_bashScriptV2.sh
[ ! -f /home/ld32/rcbio/flag/Job_bash_bashScriptV2.sh.success ] && exit 1 || exit 0

#sbatch command:
#sbatch -p short -t 10 -c 1 --mail-type=FAIL --nodes=1  -J Job_bash_bashScriptV2.sh -o /home/ld32/rcbio/flag/Job_bash_bashScriptV2.sh.out -e /home/ld32/rcbio/flag/Job_bash_bashScriptV2.sh.out /home/ld32/rcbio/flag/Job_bash_bashScriptV2.sh.sh

# Submitted batch job 43976
