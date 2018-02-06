#!/bin/sh

usage() { echo -e "Usage :\n${0##*/} [-r species_index (required if no -b. Such as: dm3, dm6, mm10, hg18, hg19 or hg38. Let us know if you need other references)] [-b BWAIndexWithPath(required if no -r, don't need this if -r is given)]"; exit 1;} 

while getopts ":r:" o; do
    case "${o}" in
        r)
            r=${OPTARG}
            ;;
        b)
            b=${OPTARG}
            ;;
    
    esac
done

module load bwa/0.7.8 picard/2.8.0 samtools/0.1.19

SHARED_DATABASES=/n/shared_db/ignome/03032016/

if [ -z "${r}" ]; then
    if [ ! -z "$b" ]; then 
        index=$b
    else
        usage
    fi
else 
  case "$r" in
    "mm10")index="$SHARED_DATABASES/Mus_musculus/UCSC/mm10/Sequence/BWAIndex/genome.fa"
    ;;
    
    "dm3") index="$SHARED_DATABASES/Drosophila_melanogaster/UCSC/dm3/Sequence/BWAIndex/genome.fa"
    ;;
    
    "dm6") index="$SHARED_DATABASES/Drosophila_melanogaster/UCSC/dm6/Sequence/BWAIndex/genome.fa"
    ;;
    
    "hg18") index="$SHARED_DATABASES/Homo_sapiens/UCSC/hg18/Sequence/BWAIndex/genome.fa"
    ;;
    
    "hg19") index="$SHARED_DATABASES/Homo_sapiens/UCSC/hg19/Sequence/BWAIndex/genome.fa"
    ;;
     
    "hg38") index="$SHARED_DATABASES/Homo_sapiens/UCSC/hg38/Sequence/BWAIndex/genome.fa"
    ;;
    
    *)  echo "Index '$r' is not supported. Please email rchelp@hms.harvard.edu for help."; exit
    ;;
    
  esac
fi

[ -d group2 ] || { echo group2 is not found. You need at least two groups to run this pipeline; exit 1; }

[ -f $index.ann ] || { echo BWA Index not find: $index.ann; exit 1; }

mkdir -p bwaOut

pwdhere=`pwd`

for group in `ls -v -d group*/|sed 's|[/]||g'`; do
    echo working on group:  $group
    cd $pwdhere/$group

    for sample in `ls -d */|sed 's|[/]||g'`; do
        echo working on sample: $sample
        inputsams=""

        mkdir -p $pwdhere/bwaOut/$group$sample
        for r1 in `ls $sample/*_1.fastq $sample/*_1.fq 2>/dev/null`; do
            readgroup=${r1#*/}
            readgroup=${readgroup%_*}
            r2=${r1%_*}_2${r1##*_1}
            echo working on readgroup: $readgroup
            [[ -f $r2 ]] || { echo -e "\n\n!!!Warning: read2 file '$r2' not exist, ignore this warning if you are working with single-end data\n\n"; r2=""; }

            #@1,0,bwa,index,sbatch -n 4 -p short -t 12:0:0 --mem 40G 
            bwa mem -M -t 4 -R "@RG\tID:$sample.$readgroupz\tPL:Illumina\tSM:$sample" $index $r1 $r2 > $pwdhere/bwaOut/$group$sample/$readgroup.sam

            inputsams="$inputsams INPUT=$pwdhere/bwaOut/$group$sample/$readgroup.sam"

        done

        #@2,1,merge
        rm $pwdhere/bwaOut/$group$sample/accepted_hits.bam 2>/dev/null;  java -Xmx12g -jar $PICARD/picard-2.8.0.jar MergeSamFiles VALIDATION_STRINGENCY=SILENT OUTPUT=$pwdhere/bwaOut/$group$sample/accepted_hits.bam  $inputsams SORT_ORDER=coordinate && samtools index $pwdhere/bwaOut/$group$sample/accepted_hits.bam

    done
    
done


