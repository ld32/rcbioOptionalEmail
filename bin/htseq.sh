#!/bin/sh

usage() { echo -e "Usage :\n${0##*/} <-p path_of_alignment, such as: bowtieOut, bwaOut or starOut, required> [-r species_index (required if no -g. Such as: dm3, dm6, mm10, hg18, hg19 or hg38. Let us know if you need other references)] [-g .gtf file with full path (required if no -r, don't need this if -r is given)]"; exit 1;}

while getopts ":p:r:g:" o; do
    case "${o}" in
        r)
            r=${OPTARG}
            ;;
        g)
            g=${OPTARG}
            ;;
        p)  
            p=${OPTARG}
            ;;
        
        *)
            usage
            ;;
    esac
done

[ -d group2 ] || { echo group2 is not found. You need at least two groups to run this pipeline; exit 1; }

# rename flag folder from earlier bowtie/star/bwa alignment run 
[ -f flag/alljobs.jid.first ] && { grep 2.1.mergeHTSeq flag/alljobs.jid.first >/dev/null || { mv -f flag flagEarlierRun; mkdir flag; } ; } 

[ -z "$p" ] && { echo Please provide the alignment folder for option -p such as -p bwaOut, -p bowtieOut or -p starOut;  usage; }

[ -d $p ] || { echo Folder not exist: $p; usage; }

module load gcc/6.2.0  python/2.7.12  htseq/0.9.1 samtools/0.1.19
SHARED_DATABASES=/n/shared_db/igenome/03032016

if [ -z "${r}" ]; then
    if [ ! -z "$g" ]; then 
        gtf=$g
    else
        usage
    fi
    
else 
  case "$r" in
     "dm3")
       gtf=" $SHARED_DATABASES/Drosophila_melanogaster/UCSC/dm3/Annotation/Genes/genes.gtf"
    ;;
    
    "dm6")
       gtf=" $SHARED_DATABASES/Drosophila_melanogaster/UCSC/dm6/Annotation/Genes/genes.gtf"
    ;;
   
    "mm10")
       gtf=" $SHARED_DATABASES/Mus_musculus/UCSC/mm10/Annotation/Genes/genes.gtf"
    ;;
        
    "hg18")
       gtf=" $SHARED_DATABASES/Homo_sapiens/UCSC/hg18/Annotation/Genes/genes.gtf"
    ;;
    
    "hg19") 
       gtf=" $SHARED_DATABASES/Homo_sapiens/UCSC/hg19/Annotation/Genes/genes.gtf"
    ;;
    
    "hg38") gtf="$SHARED_DATABASES/Homo_sapiens/UCSC/hg38/Annotation/Genes/genes.gtf"
    ;;
     
    *)  echo "Index '$r' is not supported. Please email rchelp@hms.harvard.edu for help."; exit
    ;;
    
  esac
fi

[ -f $gtf ] || { echo .gtf file not exist; exit 1; }

pwdhere=`pwd`

for group in `ls -v -d group*/|sed 's|[/]||g'`; do
    echo working on group:  $group
    cd $pwdhere/$group
    for sample in `ls -d */|sed 's|[/]||g'`; do
        echo working on sample: $sample
        i=$pwdhere/$p/$group$sample/accepted_hits.bam
        [ -f $i ] || { echo bam file not find: $i; exit 1; }
        
        #@1,0,htseq,,sbatch -p short -t 12:0:0 --mem 8G
        samtools view -bf 1 $i > $i.pe.bam && samtools sort -n $i.pe.bam $i.sorted && htseq-count -s no -f bam $i.sorted.bam $gtf > $i.read.count.txt

    done

done

module load R/3.4.1

cd $pwdhere

#@2,1,mergeHTSeq,,sbatch -p short -t 2:0:0 --mem 4G
mergeHTSeqCount.r $p


