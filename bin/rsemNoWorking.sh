#!/bin/sh

usage() { echo -e "\nUsage : $0 [-q queueName(optional, default: mcore)] [-n numberOfCore(optional, default: 4)] [-W runTime(optional, default: 24:00)] [-r species_index (required if no -b. Such as: mm10 or hg19. Let us know if you need other references)]"; exit 1;} 

while getopts ":q:n:W:r:b:f:g:t:" o; do
    case "${o}" in
        q)
            q=${OPTARG}
            ;;
        n)
            n=${OPTARG}
            ;;
        W)
            W=${OPTARG}
            ;;
        r)
            r=${OPTARG}
            ;;
        *)
            usage
            ;;
    esac
done
#echo 1 queue: $q core: $n run time: $W r: $r b: $b g: $g t: $t

[ -z "${q}" ] && q=mcore
[ -z "${n}" ] && n=4 
[ -z "${W}" ] && W=48:00

module load dev/java/jdk1.7 seq/rsem/1.2.25 seq/bowtie/1.1.1  seq/igvtools/2.3.57 

xsub="bsub -q $q -n $n -W $W"

case "$r" in
    "mm10")index="$SHARED_DATABASES/rsem_indexes/mm10"
    ;;
    "hg19") index="$SHARED_DATABASES/rsem_indexes/hg19"
    ;;
    
    *)  echo "Index '$r' is not supported. Please email rchelp@hms.harvard.edu for help."; exit
    ;;
esac


[ -d group2 ] || { echo group2 is not found. You need at least two groups to run this pipeline; exit 1; }

rm flag/*.jid 2>/dev/null  # job denpendency id file for each job

mkdir -p out

# go through sample groups
#loopStart,group
for group in `ls -d group*`; do
    echo working on group:  $group
    
    COUNTER=0 
    #loopStart,sample
    for sample in `ls -d $group/*  | xargs -n 1 basename`; do 
        echo working on sample: $sample
        
        ls  $group/$sample/*_1.fastq >/dev/null 2>&1 || ls $group/$sample/*_1.fq  >/dev/null 2>&1 || { echo Read file not found for $sample! Please make sure the fastq files are named as xxx_1.fastq, xxx_2.fastq or xxx_1.fq, xxx_2.fq;  exit 1; }
        
        reads1=""; reads2=""        
        for r1 in `ls $group/$sample/*_1.fastq $group/$sample/*_1.fq 2>/dev/null | xargs -n 1 basename`; do 
            readgroup=${r1##*/}
            readgroup=${readgroup%_*}
            r2=${r1%_*}_2${r1##*_1}
            echo working on readgroup: $readgroup
			reads1="$reads1$group/$sample/$r1," 
            
            [[ -f $group/$sample/$r2 ]] && reads2="$reads2$group/$sample/$r2," || { echo -e "\n!!!Warning: read2 file '$r2' not exist, ignore this warning if you are working with single-end data\n"; }
            
        
        done
        
        #echo reads1 $reads1
        #echo reads2 $reads2         
        
        [ -z $reads2 ] || reads1="--paired-end $reads1"
                 
        #@1,0,expression:    
        rsem-calculate-expression -p $n -q  --bowtie-chunkmbs 200 --output-genome-bam ${reads1%,} ${reads2%,} $index out/$group.$sample 
        
        
        #@2,1,wigTranscript:  
        rsem-bam2wig out/$group.$sample.transcript.sorted.bam out/$group.$sample.transcript.wig out/$group.$sample.transcript && igvtools toTDF out/$group.$sample.transcript.wig out/$group.$sample.transcript.tdf $index.transcripts.genome > out/$group.$sample.igvtools.out 
        
        #@3,1,wigGenome:  
        rsem-bam2wig out/$group.$sample.genome.sorted.bam out/$group.$sample.genome.wig out/$group.$sample.genome 
        
        #@4,1,plot:  
        rsem-plot-model out/$group.$sample out/$group.$sample.plot.pdf
        
        COUNTER=$((COUNTER + 1))          
        generesults="$generesults out/$group.$sample.genes.results"
        isoresults="$isoresults out/$group.$sample.isoforms.results" 
    
    #loopEnd            
    done 
    [[ "$COUNTER" -eq 0 ]] || counts="$counts,$COUNTER"

#loopEnd    
done 

#@5,1,iso:    
rsem-generate-data-matrix $isoresults > iso.count.matrix &&        rsem-run-ebseq --ngvector $index.ngvec iso.count.matrix ${counts#,} all.iso.results && rsem-control-fdr all.iso.results 0.05 final.iso.results.fdr.txt

#@6,1,gene:
rsem-generate-data-matrix $generesults > gene.count.matrix &&  rsem-run-ebseq gene.count.matrix ${counts#,} all.gene.results && rsem-control-fdr all.gene.results 0.05 final.gene.results.fdr.txt

