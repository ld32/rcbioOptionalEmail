#!/bin/sh

usage() { echo -e "\nUsage: \n$0 <bash_script_v2.sh> <bash_script_v3.sh> <sbatch options, such as: \"sbatch -p medium -t 24:0:0 -n 4\" . Notice: it should be double quoted.> <useTmp/noTmp>"; exit 1; } 

[[ "$3" != sbatch* ]] && usage

[ -f "$1" ] || { echo bash script file not exist: $1. Exiting...; usage; }

core=${3#*-n }; core=${core%% *}; echo core: $core; [ "$core" -eq "$core" ] || { echo core is not number; exit 1; } 

run="$2"

echo converting $1 to $2


cat <<EOT > $run
    #!/bin/sh 

    echo Running \$0 \$@       

    xsub=\"$3\" 

    stamp=\$(date -d \"today\" +\"%Y%m%d%H%M\")
    
    mkdir -p flag

    if [ -f flag/alljobs.jid ]; then
        
        checkJobsSlurm  flag/alljobs.jid 

        [ \$? == 1 ] && exit 0;" 
    fi

    cwd=\`realpath ./flag\`

    #echo rm flag/*.failed flag/*.killed 2>/dev/null

    [ -f flag/alljobs.jid ] && mv flag/alljobs.jid flag/alljobs.jid.old

    printf \"%-10s   %-20s   %-10s\n\" job_id depend_on job_flag > flag/alljobs.jid
  
    echo ---------------------------------------------------------
EOT


[ "$4" == "useTmp" ] && echo ". $(dirname $0)/rcUtils.sh" >> $run

IFS=$'\n'
for t in `cat $1`; do 
    #echo original row: .$t. 
    
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
        ref=${arrIN[3]%%:*} # references files which need to rsync to tmp space
	    [ -z "$ref" ] || ref=.$ref
        
        echo "$t" >> $run
    
        #echo step is $step
        #echo de is $de
        #echo name is $name
        cmd=""
    elif [ ! -z $findjob ]; then
        [[ "$i" == \#* ]] && { findjob=""; echo "$t" >> $run; continue; }
        echo
        echo find job: 
        echo $i
        
        # command maybe multiple lines
        if [[ "$i" == *\\ ]]; then  # command have multiple lines, at end of the line is a \
          #echo check i: $i
          cmd="$cmd ${i%\\}"
          #echo cmd0: $cmd
          continue
        fi
        space="$oldspace"
        cmd="$cmd ${i%\\#*}" # remove comments
        cmd="${cmd# }"       # remove leading space
        
        #echo cmd1:$cmd
        findjob=""
        [ ! -z ${find[$step]} ] && { echo job step id should be unique: $step exiting...; exit 1; } ||  find[$step]=yes
                
        cmd=${cmd//\\/\\\\}    # escape back slash
        
        #echo cmd2 is $cmd
        cmd=${cmd//\"/\\\"} # escape double quotes
        
        #echo cmd3 is $cmd
        
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
        [ "$4" == "useTmp" ] && useTmp="reference: $ref"
	    cat <<EOT >> $run
	        ${space}echo; echo step: $step, depends on: $de, job name: $name, flag: $name$loper $useTmp
	        ${space}flag=${step}.$de.${name}${loper}
            ${space}flag=\${flag//\//_}   # replace path / to _ 
        
            ${space}deps=\"\" 
EOT
        if [[ "$de" != "0" ]]; then
            if [[ "$de" == *\.* ]]; then
                for dep in ${de//\./ }; do
                    [ -z $sameloop ] && echo "${space}deps=\"\$deps \${jobID[$dep]}\"" >> $run || echo "${space}deps=\"\$deps \${jobIDs[$dep]}\"" >> $run
                done
            else    
                [ -z $sameloop ] && echo "${space}deps=\"\${jobID[$de]}\""  >> $run || echo "${space}deps=\"\${jobIDs[$de]}\""  >> $run 
            fi
        fi 
        
                
        # escape double quota. bsub does not need this? not sure
        cmd=${cmd//\"/\\\\\"}

        # replace space with ., if the job depends on something
        echo "${space}[ -z \"\${deps// /}\" ] && deps=null || deps=\${deps// /.}" >> $run

        #echo "echo command is: bsubRun \$xsub -flag \$deps \$cwd \$flag \"$cmd\" " >> $run
        [[ "$4" == "useTmp" && ! -z "$ref" ]] && echo "${space}setPath $ref" >> $run && cmd="rsyncToTmp ${ref//./ $}; $cmd"

        cmd="{ $cmd; } && touch \$cwd/\$flag.success || touch \$cwd/\$flag.failed"      
        
        # echo exit >> $run   # only test for one job
        echo "${space}id=\$(sbatchRun \$xsub -flag \$deps \$cwd \$flag \"srun -n 1 bash -c \\\"$cmd\\\"\")"   >> $run
       
        [[ "$4" == "useTmp" && ! -z "$ref" ]] && echo "${space}setPathBack $ref" >> $run 
        cat <<EOT >> $run
            #echo id is: \$id 
                    
            ${space}if [ -z \"\$id\" ]; then
            ${space}    echo  job \$flag is not submitted
            ${space}    jobID[$step]=\"\"
            ${space}else
            #${space}    touch \$cwd/\$flag.submitted 
        
            ${space}    alljobs=\"\$alljobs \$id\"
            ${space}    printf \"%-10s  %-20s  %-10s\n\" \$id \$deps \$flag >> \$cwd/alljobs.jid
EOT
        
        # tell this is out of the loop for the depending job (de), so that we clear the job id list for the next step with depends on 'de'
        if [[ "$de" == *\.* ]]; then
                for dep in ${de//\./ }; do
                    [ -z $sameloop ] && echo "${space}    startNewLoop[$dep]=\"no\""  >> $run || echo "    ${space}startNewLoop[$dep]=\"\""  >> $run 
                done
        else    
                [ -z $sameloop ] && echo "${space}    startNewLoop[$de]=\"no\""  >> $run || echo "    ${space}startNewLoop[$de]=\"\""  >> $run 
        fi
        cat <<EOT >> $run
            ${space}    [ -z \${startNewLoop[$step]} ] && jobIDs[$step]=\"\" && startNewLoop[$step]=\"no\" 
            ${space}    jobID[$step]=\$id"  >> $run
            ${space}    jobIDs[$step]=\${jobIDs[$step]}.\$id
        
            ${space}fi 
EOT
        IFS=''
        cmd=""
          
    else
        echo "$t" >> $run
        #findjob=""
        if [[ "$i" == \#loopEnd* ]]; then
            echo find loopend: $i
            #space=${space%    }
            loper=${loper%.\$*}
        elif [[ "$i" == \#loopStart* ]]; then
            echo; echo find loopStart: $i
            a=`echo $i | xargs`  # remove leading space and tailing space
            a=${a#*,}
            #echo a is .$a.
            loper="$loper.\$$a"
            #echo new loper: $loper
        fi    
    fi   
done
[ $? == 1 ] && exit 0;

# go back to the original folder

cat <<EOT >> $run

    cd \$cwd/..

    echo all submitted jobs: 
    cat flag/alljobs.jid 
    echo ---------------------------------------------------------
    [ -f flag/alljobs.jid.first ] || cp flag/alljobs.jid flag/alljobs.jid.first 
EOT

chmod a+x $run 

echo all done

