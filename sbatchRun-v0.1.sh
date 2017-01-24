#!/bin/sh

#set -x 
# to call this:                                        1                 2                                3         4          5        
# sbatchRun sbatch -p priority -n 1 -t 2:0:0 -flag   <depend_job_ids>    <absolute_path_to_flag_folder>   <my_job>  <"my_application para1   para2"> 

# This script will submit a sbatch job named "my_job", run commands in "my_job.sh", output to "my_job.out".
# When the job runs, it creates "my_job.start", "my_job.failed" if fails, "my_job.success" if successes. 
# If this job fails, it kills all jobs which depends on this job.

# If this is a re-run it also delete the downstream jobs' .success flag file, so that user is not asked to if he want to re-run these downstream jobs 
para=( "$@" )

for i in "${!para[@]}"; do 
    #echo -e "$i ${para[$i]}" >&2
    if [ "${para[$i]}" == "-flag" ]; then
        b="${para[@]:0:$i}"
        d=${para[$i+1]}
        w=${para[$i+2]}
        f=${para[$i+3]}
        c="${para[$i+4]}"
        # check the application is available (module is loaded)
        #type "${para[$i+2]}"  >/dev/null 2>&1 || { echo >&2 "Error: command '${para[$i+2]}' could not be found, did you forget to load the module? Aborting."; exit 1; }
    fi
done


#echo b: $b  >&2
#echo d: $d   >&2
#echo w: $w  >&2
#echo f: $f  >&2
#echo c: $c >&2

# escape back slash
c=${c/\\/\\\\}
#echo c: $c >&2

d=${d#.}

#echo d1 ${d/\./} >&2
if [[ $d == null ]]; then
    deps=""
    echo depend on no job >&2
elif [[ $d == ${d/\./} ]]; then
    echo depend on single job >&2
    deps="--dependency=afterok:${d/\./}"
else
    echo depend on multiple jobs >&2
    for t in ${d//\./ }; do
       deps="$deps:$t"
    done
    [ ! -z "$deps" ] && deps="--dependency=afterok$deps"
fi

#echo checking file: $w/$f.success  >&2

#touch $w/$f.success

if [ -f $w/$f.success ]; then
    read -p "$f was done before, do you want to re-run for $sample (y)?:" x </dev/tty
	if [[ "$x" == "y" ]]; then
        rm $w/$f.success 
        
        # remove downstream jobs' .success flag
        findDownStreamJobs.sh $w $f removeSuccessFlag >&2 
        echo command is: findDownStreamJobs.sh $w $f removeSuccessFlag >&2 
        echo Will re-run the down stream steps even if they are done before. >&2 
    else
        exit; 
    fi  
fi

 
rm -f $w/$f.start $w/$f.failed $w/$f.killed  2>/dev/null 

echo -e "#!/bin/bash\n$other\n#Commands:\ntouch $w/$f.start\n$c && touch $w/$f.success || { touch $w/$f.failed; exit 1; }" > $w/$f.sh  

# if this job fails, kill all the jobs which depends on this one
# remove the line if you don't need to kill jobs automatically
# don't need to kill for slurm, the jobs will be killed automatically
#dep="-Ep \"${0%/*}/findDownStreamJobs.sh $w $f \""; 

cmd="$b  --kill-on-invalid-dep=yes --nodes=1  $dep --mem-per-cpu=10G $deps -J $f --mail-type=ALL -o $w/$f.out -e $w/$f.out --wrap=\"sh $w/$f.sh\"" 

#cmd="$b $deps -J $f -B -N $dep sh $w/$f.sh" 

[[ "$deps" != *&&* ]] && output="$($cmd)" || output="$(eval $cmd)"
#echo bsub output: $output >&2
#echo $output | head -n1 | cut -d'<' -f2 | cut -d'>' -f1 
echo ${output##* }
#  tee -a $f.jid; #here don't need to write the id any more.

#cat $f.sh >&2
echo -e "# $cmd \n# $output\n" >&2

echo -e "\n# $output" >> $w/$f.sh
