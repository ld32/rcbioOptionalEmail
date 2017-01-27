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
            bwa mem -M -t 4 -R "@RG\tID:$sample.$readgroupz\tPL:Illumina\tSM:$sample" $index $r1 $r2 > $sample/$readgroup.sam 
            
            inputsams="$inputsams INPUT=$sample/$readgroup.sam"
   
        #loopEnd 	
		done
        
        #@2,1,merge:merge bam files 
        java -Xmx12g -jar $PICARD/picard-2.8.0.jar MergeSamFiles VALIDATION_STRINGENCY=SILENT OUTPUT=$sample.bam  $inputsams SORT_ORDER=coordinate && samtools index $sample.bam
        
        bams="$bams -I $group/$sample.bam"
      
    #loopEnd 
    done 
#loopEnd     
done

cd $pwdhere

echo $bams > bam.list
