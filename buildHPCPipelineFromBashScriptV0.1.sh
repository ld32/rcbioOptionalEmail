#!/bin/sh

usage() { echo -e "\nUsage: \n$0 <old_bash_script.sh> <new_hpc_pipeline.sh> <bsub options, such as: \"bsub -W 30:00 -q mcore -n 4\" . Notice: it should be double quoted.>"; exit 1; } 

[[ "$3" != sbatch* &&  "$3" != bsub* ]] && usage

[ -f "$1" ] || { echo bash script file not exist: $1. Exiting...; usage; }

run="$2"

echo converting $1 to $2

echo "#!/bin/sh" > $run  

echo "echo Running \$0 \$@"  >> $run        

echo "xsub=\"$3\"" >> $run

echo "stamp=\$(date -d \"today\" +\"%Y%m%d%H%M\")" >> $run
    
echo mkdir -p flag >> $run

echo "if [ -f flag/alljobs.jid ]; then" >> $run
[[ "$3" == sbatch* ]] && echo "    checkJobsSlurm  flag/alljobs.jid " >> $run
[[ "$3" == bsub* ]] &&   echo "    checkJobs  flag/alljobs.jid " >> $run



echo "    [ \$? == 1 ] && exit 0;" >> $run
echo "fi" >> $run


echo "cwd=\`realpath ./flag\`" >> $run

#echo "rm flag/*.failed flag/*.killed 2>/dev/null" >> $run

echo "[ -f flag/alljobs.jid ] && mv flag/alljobs.jid flag/alljobs.jid.old"  >> $run 

echo "printf \"%-10s   %-20s   %-10s\n\" job_id depend_on job_flag > flag/alljobs.jid" >> $run 
echo "echo ---------------------------------------------------------" >> $run
#echo "for i in \$alljobs; do" >> $run
#echo "    printf \"%-10s  %-20s  %-10s\n\" \$i \${mmap[\$i]} \${nmap[\$i]} >> flag/alljobs.jid"  >> $run 
#echo done >> $run

#echo "declare -A nmap" >> $run
#echo "declare -A mmap" >> $run

IFS=$'\n'
for t in `cat $1`; do 
    #echo originalis:.$t. 
    
    # get space and the real command
    space=""
    for (( i=0; i<${#t}; i++ )); do
        chr=${t:$i:1}
        if [[ "$chr" == " " || "$chr" == $'\t' ]]; then
            space="$space$chr"
        else
           # i is the real command
           i="${t:$i:1000}"
           break
        fi
    done
    
    
    # this maybe a new command, remember the space
    [ -z "$cmd" ] && oldspace="$space"
    
    if [[ "$i" == \#@* ]]; then
        echo
        echo find job marker: 
        echo $i
        
        findjob="yes"
        x=${i#*\#@}
        IFS=',' read -a arrIN <<< "$x"
        IFS=''
    
        step=${arrIN[0]}
       
        
        de=${arrIN[1]}   
        name=${arrIN[2]%%:*}
        echo "$t" >> $run
    
        #echo step is $step
        #echo de is $de
        #echo name is $name
        cmd=""
    elif [ ! -z $findjob ]; then
        #findjob=""
        [[ "$i" == \#* ]] && { findjob=""; echo "$t" >> $run; continue; }
        echo
        echo find job: 
        echo $i
        #make sure the job id is unique in the pipeline
        #[ ! -z ${find[$step]} ] && { echo job step id should be unique: $step exiting...; exit 1; } ||  find[$step]=yes
        # command maybe multiple lines
        if [[ "$i" == *\\ ]]; then
          #echo check i: $i
          cmd="$cmd ${i%\\}"
          echo cmd0: $cmd
          continue
        fi
        space="$oldspace"
        cmd="$cmd ${i%\\#*}" # remove comments
        cmd="${cmd# }"   # remove leading space
        #echo cmd1:$cmd
        findjob=""
        [ ! -z ${find[$step]} ] && { echo job step id should be unique: $step exiting...; exit 1; } ||  find[$step]=yes
                
        cmd=${cmd//\\/\\\\}    # escape back slash
        
        #echo cmd is $cmd
        cmd=${cmd//\"/\\\"} # escape double quotes
        
        #echo cmd is$cmd
        
        cloper[$step]=$loper
        sameloop=no
        IFS=' '
        if [[ "$de" == "0" ]]; then
            sameloop=""
        elif [[ "$de" == *\.* ]]; then
            
            for dep in ${de//\./ }; do
                if [[ ${cloper[$step]} == ${cloper[$dep]} ]]; then
                    sameloop=""
                else
                    sameloop="no"
                    break
                fi
            done
        elif [[ ${cloper[$step]} == ${cloper[$de]} ]]; then
            sameloop=""
        fi
                      
        #echo "debug: ${cloper[$step]} ${cloper[$de]} .$sameloop. " >> $run
        echo -e "\n${space}echo; echo step: $step, depends on: $de, job name: $name, flag: $name$loper  " >> $run
        echo "${space}flag=${step}.$de.${name}${loper}"  >> $run
        echo "${space}flag=\${flag//\//_}" >> $run   # replace path / to _ 
        
        echo -e "${space}deps=\"\""  >> $run
        if [[ "$de" != "0" ]]; then
            if [[ "$de" == *\.* ]]; then
                for dep in ${de//\./ }; do
                    [ -z $sameloop ] && echo "${space}deps=\"\$deps \${jobID[$dep]}\"" >> $run || echo "${space}deps=\"\$deps \${jobIDs[$dep]}\"" >> $run
                    
                    #echo  "${space}findDownStreamJobs.sh \$cwd $dep removeSuccessFlag" >> $run 
                done
            else    
                [ -z $sameloop ] && echo "${space}deps=\"\${jobID[$de]}\""  >> $run || echo "${space}deps=\"\${jobIDs[$de]}\""  >> $run 
                
                #echo  "${space}findDownStreamJobs.sh \$cwd $de removeSuccessFlag" >> $run 
            fi
            
            # check if any of the upstream job is being rerun, if yes, delete the .success file for this job
        #    echo "${space}tempids=\"\"" >> $run           
        #    if [ -z $sameloop ]; then
        #        if [[ "$de" == *\.* ]]; then
        #            for dep in ${de//\./ }; do
        #                echo "${space}tempids=\"\$tempids \${depes[$dep]}\"" >> $run    
        #            done
        #        else
        #            echo "${space}tempids=\"\$tempids \${depes[$de]}\"" >> $run                
        #        fi
        #    else
        #        if [[ "$de" == *\.* ]]; then
        #            for dep in ${de//\./ }; do
        #                echo "${space}tempids=\"\$tempids \${depesm[$dep]}\"" >> $run     
        #            done
        #        else
        #            echo "${space}tempids=\"\$tempids \${depesm[$de]}\"" >> $run             
        #        fi
        #    fi
        #    echo "${space}sorted_unique_ids=\$(echo \"\$tempids\" | tr ' ' '\n' | sort -u | tr '\n' ' ') " >> $run
        #    echo "${space}for d in \"\"\$sorted_unique_ids\"\"; do " >> $run 
        #    echo "echo checking file \$cwd/\${d%.*}.redo" >> $run
        #    #echo "[[ -f \$cwd/\${d%.*}.redo && -f \$cwd/\$flag.success ]] && echo removed: \$cwd/\$flag.success " >> $run
        #    echo "${space}    [[ -f \$cwd/\${d%.*}.redo && -f \$cwd/\$flag.success ]] && rm \$cwd/\$flag.success " >> $run
        #    
        #    echo "${space}done " >> $run     
        fi 
  
       
        
        #echo "${space}for i in \${deps//./ }; do" >> $run
        #echo "${space}    findUpStreamJobs.sh \$cwd \$i \$flag" >> $run
        #echo "${space}done" >> $run         
        
        # replace space with ., if the job depends on something, remove .success flag if exist
        echo "${space}[ -z \"\${deps// /}\" ] && deps=null || deps=\${deps// /.}" >> $run

        #echo "echo command is: bsubRun \$xsub -flag \$deps \$cwd \$flag \"$cmd\" " >> $run
        
        [[ "$3" == sbatch* ]] && echo "${space}id=\$(sbatchRun \$xsub -flag \$deps \$cwd \$flag \"$cmd\")"   >> $run
        [[ "$3" == bsub* ]] && echo "${space}id=\$(bsubRun \$xsub -flag \$deps \$cwd \$flag \"$cmd\")"   >> $run

        #echo "echo id is: \$id ">> $run
        
        echo "${space}if [ -z \"\$id\" ]; then"  >> $run
        echo "${space}    echo  job \$flag is not submitted"  >> $run
        echo "${space}    jobID[$step]=\"\"" >> $run 
        echo "${space}else"  >> $run
        #echo "${space}    touch \$cwd/\$flag.submitted" >> $run 
        
        echo "${space}    alljobs=\"\$alljobs \$id\"" >> $run 
        #echo "${space}    nmap[\$id]=\$flag"  >> $run        
        #echo "${space}    mmap[\$id]=\$deps"  >> $run
        echo "${space}    printf \"%-10s  %-20s  %-10s\n\" \$id \$deps \$flag >> \$cwd/alljobs.jid"  >> $run 
        
        # tell this is out of the loop for the depending job (de), so that we clear the job id list for the next step with depends on 'de'
        if [[ "$de" == *\.* ]]; then
                for dep in ${de//\./ }; do
                    [ -z $sameloop ] && echo "${space}    startNewLoop[$dep]=\"no\""  >> $run || echo "    ${space}startNewLoop[$dep]=\"\""  >> $run 
                done
        else    
                [ -z $sameloop ] && echo "${space}    startNewLoop[$de]=\"no\""  >> $run || echo "    ${space}startNewLoop[$de]=\"\""  >> $run 
        fi
        
        #echo "${space}    [ -z \${startNewLoop[$step]} ] && jobIDs[$step]=\"\" && startNewLoop[$step]=\"no\" && echo starting new loop for $step!; echo "  >> $run
 
        echo "${space}    [ -z \${startNewLoop[$step]} ] && jobIDs[$step]=\"\" && startNewLoop[$step]=\"no\" " >> $run 
         
         
        echo "${space}    jobID[$step]=\$id"  >> $run
        
       
        
        echo "${space}    jobIDs[$step]=\${jobIDs[$step]}.\$id"  >> $run
        
        
        
        
        #echo "${space}    depes[$step]=\$flag.jid" >> $run 
        #echo "${space}    depesm[$step]=\"\${depesm[$step]} \$flag.jid\"" >> $run 
        #if [[ "$de" == *\.* ]]; then             
        #    for dep in ${de/\./ }; do
        #       echo "${space}    depes[$step]=\"\${depes[$step]} \${depes[$dep]}\"" >> $run  
        #        echo "${space}    depesm[$step]=\"\${depesm[$step]} \${depesm[$dep]}\"" >> $run
        #    done
        #else
        #    echo "${space}    depes[$step]=\"\${depes[$step]} \${depes[$de]}\"" >> $run 
        #    echo "${space}    depesm[$step]=\"\${depesm[$step]} \${depesm[$de]}\"" >> $run             
        #fi
        #if [ -z $sameloop ]; then
        #    echo "${space}    sorted_unique_ids=\$(echo \"\${depes[$step]}\" | tr ' ' '\n' | sort -u | tr '\n' ' ') " >> $run
        #else
        #    echo "${space}    sorted_unique_ids=\$(echo \"\${depesm[$step]}\" | tr ' ' '\n' | sort -u | tr '\n' ' ') " >> $run
        #fi
        #echo "${space}    for d in \"\"\$sorted_unique_ids\"\"; do " >> $run 
        #echo "${space}       echo \$id \${mmap[\$id]} \${nmap[\$id]} >> \$cwd/\$d" >> $run  
        #echo "echo writing id: \$id to \$cwd/\$d " >> $run
        #echo "${space}    done " >> $run  
        echo "${space}fi" >> $run 
        IFS=''
        cmd=""
             
    else
        echo "$t" >> $run
        #findjob=""
        if [[ "$i" == \#loopEnd* ]]; then
            echo
            echo find loopend: $i
            #space=${space%    }
            loper=${loper%.\$*}
        #fi
        elif [[ "$i" == \#loopStart* ]]; then
            echo
            echo find loopStart: $i
            #loopStart='find'
            #space="${space}    "
            a=${i#*,}
            loper="$loper.\$$a"
            
            
        fi    
    fi   
done
[ $? == 1 ] && exit 0;

# go back to the original folder
echo "cd \$cwd/.." >> $run

echo "echo all submitted jobs: " >> $run
#echo "printf \"%-10s   %-20s   %-10s\n\" job_id depend_on job_flag" >> $run 
#echo "echo ---------------------------------------------------------" >> $run
#echo "for i in \$alljobs; do" >> $run
#echo "    printf \"%-10s | %-20s | %-10s\n\" \$i \${mmap[\$i]} \${nmap[\$i]}" >> $run 
#echo done >> $run
echo "cat flag/alljobs.jid" >> $run 

#echo "printf \"%-10s   %-20s   %-10s\n\" job_id depend_on job_flag > flag/alljobs.jid" >> $run 
echo "echo ---------------------------------------------------------" >> $run
#echo "for i in \$alljobs; do" >> $run
#echo "    printf \"%-10s  %-20s  %-10s\n\" \$i \${mmap[\$i]} \${nmap[\$i]} >> flag/alljobs.jid"  >> $run 
#echo done >> $run
#echo sleep 5 >> $run

#echo echo bjobs -w output: >> $run
#echo bjobs -w >> $run

chmod a+x $run 

echo "[ -f flag/alljobs.jid.first ] || cp flag/alljobs.jid flag/alljobs.jid.first " >> $run 

echo all done

