#!/bin/sh


usage() { echo -e "Usage :\n${0##*/} [-r species_index (required if no -b. Such as: dm3, dm6, mm10, hg18, hg19 or hg38. Let us know if you need other references)] [-b bowtie2IndexWithPath(required if no -r, don't need this if -r is given)]"; exit 1;} 

while getopts ":r:b:" o; do
    case "${o}" in
        r)
            reference=${OPTARG}
            ;;
        b)
            bowtieIndex=${OPTARG}
            ;;
    
    esac
done

module load bowtie2/2.2.9 samtools/0.1.19 
SHARED_DATABASES=/n/shared_db/igenome/03032016/

if [ -z "${reference}" ]; then
    if [ ! -z "$bowtieIndex" ]; then 
        bowtie2-inspect -n ${bowtieIndex} &> /dev/null || { echo -e "Error: \ngenome bowtie2 index could not be found: $bowtieIndex"; usage; } 
        index=$bowtieIndex
    else
        usage
    fi
    
else 
  case "$reference" in
    "mm10")index="$SHARED_DATABASES/Mus_musculus/UCSC/mm10/Sequence/Bowtie2Index/genome"
    ;;
    
    "dm3") index="$SHARED_DATABASES/Drosophila_melanogaster/UCSC/dm3/Sequence/Bowtie2Index/genome"
    ;;
    
    "dm6") index="$SHARED_DATABASES/Drosophila_melanogaster/UCSC/dm6/Sequence/Bowtie2Index/genome"
    ;;
    
    "hg18") index="$SHARED_DATABASES/Homo_sapiens/UCSC/hg18/Sequence/Bowtie2Index/genome"
    ;;
    
    "hg19") index="$SHARED_DATABASES/Homo_sapiens/UCSC/hg19/Sequence/Bowtie2Index/genome"
    ;;
     
    "hg38") index="$SHARED_DATABASES/Homo_sapiens/UCSC/hg38/Sequence/Bowtie2Index/genome"
    ;;
    
    *)  echo "Index '$r' is not supported. Please email rchelp@hms.harvard.edu for help."; exit
    ;;
    
  esac
fi

[ -d group2 ] || { echo group2 is not found. You need at least two groups to run this pipeline; exit 1; }

mkdir -p bowtieOut

for group in `ls -v -d group*/|sed 's|[/]||g'`; do
    echo working on group:  $group
    #loopStart,sample
    for sample in `ls -d $group/*/ | xargs -n 1 basename`; do
         
        echo working on sample: $sample
        #inputsams=""
        ls  $group/$sample/*_1.fastq* 2>/dev/null || ls $group/$sample/*_1.fq*  2>/dev/null || { echo Read file not found for $sample! Please make sure the fastq files are named as xxx_1.fastq, xxx_2.fastq or xxx_1.fq, xxx_2.fq;  exit 1; }
        read1=""
        read2=""
        
        for r1 in `ls $group/$sample/*_1.fastq* $group/$sample/*_1.fq* 2>/dev/null | xargs -n 1 basename`; do 
            #echo r1 is $r1
            readgroup=${r1%_*}
            echo working on readgroup: $readgroup
            r2=${r1%_*}_2${r1##*_1}
            #echo R2 is $r2
            
			[[ -f $group/$sample/$r2 ]] && r2="$group/$sample/$r2" || { echo -e "\n\n!!!Warning: read2 file '$r2' not exist, ignore this warning if you are working with single-end data\n\n"; r2=""; }
            read1="$read1,$group/$sample/$r1"; [ -z $r2 ] || read2="$read2,$r2"
        done
        
        #@1,0,bowtie2,index,sbatch -n 4 -p short -t 12:0:0 --mem 40G 
        rm -r bowtieOut/$group$sample 2>/dev/null ; mkdir -p bowtieOut/$group$sample;  bowtie2 -p 4 -x $index -1 ${read1#,} -2 ${read2#,} | samtools view -bS - > bowtieOut/$group$sample/accepted_hits.bam && samtools sort bowtieOut/$group$sample/accepted_hits.bam bowtieOut/$group$sample/accepted_hits_sorted &&  samtools index bowtieOut/$group$sample/accepted_hits_sorted.bam

    done 
done

