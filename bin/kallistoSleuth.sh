#!/bin/sh
usage() { echo -e "\nUsage : `basename $0` <-i transcriptome index (required: currently we have: GRCh38, GRCm38, BDGP6 and GRCz10)> [ -l fragmentLength, optional. Only need this for single end reads]\n"; exit 1;} 

while getopts ":i:l:" o; do
    case "${o}" in
        i)
            i=${OPTARG}
            ;;

        l)  l=${OPTARG}
	    ;;
        *)
            usage
            ;;
    esac
done

[ -z "$i" ] && usage

module load gcc/6.2.0 kallisto/0.43.1  R/3.4.1 

SHARED_DATABASES=/n/shared_db

# set the correct index and gene annotation file
if [ -f "${i}" ]; then
   index=$i
else 
  case "$i" in
    "GRCh38")index="$SHARED_DATABASES/GRCh38/91/kallisto/0.43.1/GRCh38"
            s=hsapiens_gene_ensembl
    ;;
    
    "GRCm38")index="$SHARED_DATABASES/GRCm38/91/kallisto/0.43.1/GRCm38"
            s=mmusculus_gene_ensembl
    ;;
    
    "BDGP6")index="$SHARED_DATABASES/BDGP6/91/kallisto/0.43.1/BDGP6"
            s=dmelanogaster_gene_ensembl
    ;;
    
    "GRCz10")index="$SHARED_DATABASES/GRCz10/91/kallisto/0.43.1/GRCz10"
            s=drerio_gene_ensembl
    ;;
  
  
    *)  echo "Index '$i' is not supported. Please email rchelp@hms.harvard.edu for help."; usage
    ;;
  esac
fi

[ -d group2/ ] || { echo group2 is not found. You need at least two groups to run this pipeline; exit 1; }

echo -e "sample\tgroup\tpath" > sample.lst 

mkdir -p kallistoOut

#loopStart,group
for group in `ls -d group*/|sed 's|[/]||g'`; do
  
    echo working on group:  $group
    #loopStart,sample
    for sample in `ls -d $group/*/ | xargs -n 1 basename`; do
         
        echo working on sample: $sample
        fileList=`ls  $group/$sample/*{_1.fastq,_1.fq,_r1.fastq,_r1.fq,_R1.fastq,_R1.fq}* 2>/dev/null`
        echo fileList: $fileList
        
        [ -z "$fileList" ] && { echo Read file not found for $group/$sample! Please make sure the fastq files are named as xxx_1.fastq, xxx_2.fastq or xxx_1.fq, xxx_2.fq or xxx_r1.fastq, xxx_r2.fastq or xxx_R1.fq, xxx_R2.fq or xxx_1.fastq.gz, xxx_2.fastq.gz or xxx_1.fq.gz, xxx_2.fq.gz or xxx_r1.fastq.gz, xxx_r2.fastq.gz or xxx_R1.fq.gz, xxx_R2.fq.gz;  exit 1; }
        
        reads=""
        for r1 in `echo $fileList | xargs -n 1 basename`; do 
            #echo working on read one file: $r1
            
            readgroup=${r1%_*}
            echo working on readgroup: $readgroup
            
            if [  -z "$l" ] ; then
                ext=${r1##*_}; ext=${ext/1/2}
                r2=${readgroup}_$ext
                reads="$reads  $group/$sample/$r1 $group/$sample/$r2"
            else                       
                reads="$reads  $group/$sample/$r1"     
            fi
        done
        
        echo -e "$sample\t$group\tkallistoOut/$group$sample" >> sample.lst 
        
        [ -z "$l" ] ||  reads="--single $reads -l $l"
        
        #@1,0,kallisto2,,sbatch -p short --mem 32G -n 4 -t 2:0:0 
        rm -r kallistoOut/$group$sample 2>/dev/null ; kallisto quant -o kallistoOut/$group$sample -b 100 -t 4 -i $index $reads
          
        
    #loopEnd 
    done 
#loopEnd     
done

#@2,1,sleuth 
Rscript sleuth.r $s
