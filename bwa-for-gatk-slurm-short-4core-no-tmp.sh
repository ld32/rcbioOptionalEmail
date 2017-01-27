#!/bin/sh
echo Running $0 $@
xsub="sbatch -p short -t 2:0:0 -n 4"
stamp=$(date -d "today" +"%Y%m%d%H%M")
mkdir -p flag
if [ -f flag/alljobs.jid ]; then
    checkJobs flag/alljobs.jid 
    [ $? == 1 ] && exit 0;
fi
cwd=`realpath ./flag`
[ -f flag/alljobs.jid ] && mv flag/alljobs.jid flag/alljobs.jid.old
printf "%-10s   %-20s   %-10s\n" job_id depend_on job_flag > flag/alljobs.jid
echo ---------------------------------------------------------
#!/bin/sh
module load gcc/6.2.0 bwa/0.7.15 picard/2.8.0 samtools/0.1.19
#[ -d group2 ] || { echo group2 is not found. You need at least two groups to run this pipeline; exit 1; }
index="/n/groups/shared_databases/gatk/ftp.broadinstitute.org/bundle/2.8/b37/human_g1k_v37.fasta"
ls group* &> /dev/null || { echo sample group folder not found. Exiting...; exit; }
pwdhere=`pwd`
# go through sample groups
#loopStart,group
for group in `ls -d group*/|sed 's|[/]||g'`; do
    echo working on group:  $group
    cd $pwdhere/$group
    #loopStart,sample
    for sample in `ls -d */|sed 's|[/]||g'`; do
        echo working on sample: $sample
        inputsams=""
        #loopStart,r1
        for r1 in `ls $sample/*_1.fastq $sample/*_1.fq 2>/dev/null`; do 
            readgroup=${r1#*/}
            readgroup=${readgroup%_*}
            r2=${r1%_*}_2${r1##*_1}
            echo working on readgroup: $readgroup
			[[ -f $r2 ]] || { echo -e "\n\n!!!Warning: read2 file '$r2' not exist, ignore this warning if you are working with single-end data\n\n"; r2=""; }
            
            #@1,0,bwa:bwa command 

            echo; echo step: 1, depends on: 0, job name: bwa, flag: bwa.$group.$sample.$r1  
            flag=1.0.bwa.$group.$sample.$r1
            flag=${flag//\//_}
            deps=""
            [ -z "${deps// /}" ] && deps=null || deps=${deps// /.}
            id=$(sbatchRun $xsub -flag $deps $cwd $flag " bwa mem -M -t 4 -R \"@RG\\tID:$sample.$readgroupz\\tPL:Illumina\\tSM:$sample\" $index $r1 $r2 > $sample/$readgroup.sam ")
            if [ -z "$id" ]; then
                echo  job $flag is not submitted
                jobID[1]=""
            else
                alljobs="$alljobs $id"
                printf "%-10s  %-20s  %-10s\n" $id $deps $flag >> $cwd/alljobs.jid
                startNewLoop[0]="no"
                [ -z ${startNewLoop[1]} ] && jobIDs[1]="" && startNewLoop[1]="no" 
                jobID[1]=$id
                jobIDs[1]=${jobIDs[1]}.$id
            fi
            
            inputsams="$inputsams INPUT=$sample/$readgroup.sam"
   
        #loopEnd 	
		done
        
        #@2,1,merge:merge bam files 

        echo; echo step: 2, depends on: 1, job name: merge, flag: merge.$group.$sample  
        flag=2.1.merge.$group.$sample
        flag=${flag//\//_}
        deps=""
        deps="${jobIDs[1]}"
        [ -z "${deps// /}" ] && deps=null || deps=${deps// /.}
        id=$(sbatchRun $xsub -flag $deps $cwd $flag " java -Xmx12g -jar $PICARD/picard-2.8.0.jar MergeSamFiles VALIDATION_STRINGENCY=SILENT OUTPUT=$sample.bam  $inputsams SORT_ORDER=coordinate && samtools index $sample.bam")
        if [ -z "$id" ]; then
            echo  job $flag is not submitted
            jobID[2]=""
        else
            alljobs="$alljobs $id"
            printf "%-10s  %-20s  %-10s\n" $id $deps $flag >> $cwd/alljobs.jid
            startNewLoop[1]=""
            [ -z ${startNewLoop[2]} ] && jobIDs[2]="" && startNewLoop[2]="no" 
            jobID[2]=$id
            jobIDs[2]=${jobIDs[2]}.$id
        fi
        
        bams="$bams -I $group/$sample.bam"
      
    #loopEnd 
    done 
#loopEnd     
done
cd $pwdhere
echo $bams > bam.list
cd $cwd/..
echo all submitted jobs: 
cat flag/alljobs.jid
echo ---------------------------------------------------------
echo bjobs -w output:
bjobs -w
[ -f flag/alljobs.jid.first ] || cp flag/alljobs.jid flag/alljobs.jid.first 
