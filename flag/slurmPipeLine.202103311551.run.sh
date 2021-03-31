#!/bin/sh
echo Running $0 $@
echo module list:
module list
xsubd="sbatch -p short -t 10 -c 1"
mkdir -p flag
if [ -f flag/alljobs.jid ]; then
    checkJobsSlurm  flag/alljobs.jid 
    [ $? == 1 ] && exit 0;
fi
cwd=`realpath ./flag`
[ -f flag/alljobs.jid ] && mv flag/alljobs.jid flag/alljobs.jid.old
flag=Job_bash_bashScriptV2.sh
id=$(sbatchRun $xsubd -flag null $cwd $flag "srun -n 1 bash -e -c \"{ set -e; bash bashScriptV2.sh; } && touch $cwd/$flag.success || touch $cwd/$flag.failed\"")
printf "%-10s   %-20s   %-10s\n" job_id depend_on job_flag > flag/alljobs.jid
if [ -z "$id" ]; then
    echo Job $flag is not submitted
else
    printf "%-10s  %-20s  %-10s\n" $id null $flag >> $cwd/alljobs.jid
    echo 
fi
cp flag/alljobs.jid flag/alljobs.jid.first
