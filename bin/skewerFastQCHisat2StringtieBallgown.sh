#!/bin/sh

usage() { echo -e "\nUsage : $0 <-r species_index -a <adapter.fa> (required. Such as: dm6, mm10, hg19 or GRCz10. Let us know if you need other references)> "; exit 1;} 

while getopts ":a:r:" o; do
    case "${o}" in
        r)
            r=${OPTARG}
            ;;
        a)
            a=${OPTARG}
            ;;    
        *)
            usage
            ;;
    esac
done

[ -z $r ] && { echo Genome referenc is needed; usage; }
[ -f $a ] || { echo Adapter sequence file is missing or not exist $a; usage; }
 
module load gcc/6.2.0 skewer/0.2.2 fastqc/0.11.5 hisat2/2.1.0 stringtie/1.3.3b samtools/1.3.1 R/3.4.1  

#python/2.7.12 #htseq/0.9.1 
#R/3.4.1 

SHARED_DATABASES=/n/groups/shared_databases/
 
case "$r" in
    "dm6" )index=/n/shared_db/dm6/hisat2/2.1.0/genome
	   splice="--known-splicesite-infile /n/shared_db/dm6/hisat2/2.1.0/splicesites.txt"
	   gtf=/n/shared_db/dm6/hisat2/2.1.0/genes.gtf
    ;;	
    "mm10")index="$SHARED_DATABASES/hisat2_indexes/mm10/genome"
           splice="--known-splicesite-infile $SHARED_DATABASES/hisat2_indexes/mm10/splicesites.txt"
            gtf=" $SHARED_DATABASES/hisat2_indexes/mm10/genes.gtf"
    ;;
    "hg19") index="$SHARED_DATABASES/hisat2_indexes/hg19/genome"
            splice="--known-splicesite-infile $SHARED_DATABASES/hisat2_indexes/hg19/splicesites.txt"
            gtf=" $SHARED_DATABASES/hisat2_indexes/hg19/genes.gtf"
    ;;
    
    "GRCz10") index="$SHARED_DATABASES/hisat2_indexes/GRCz10/genome"
            splice="--known-splicesite-infile $SHARED_DATABASES/hisat2_indexes/GRCz10/splicesites.txt"
            gtf=" $SHARED_DATABASES/hisat2_indexes/GRCz10/genes.gtf"
    ;;
    
    *)  echo "Index '$r' is not supported. Please email rchelp@hms.harvard.edu for help."; usage
    ;;
esac

 
[ -d group2 ] || { echo group2 is not found. You need at least two groups to run this pipeline; exit 1; }

pwd=`pwd`

for group in `ls -v -d group*/|sed 's|[/]||g'`; do
    echo working on group:  $group
    #loopStart,sample
    for sample in `ls -d $group/*/ | xargs -n 1 basename`; do
        groups="$groups,${group#group}"; 
        echo working on sample: $sample, read1 files: 
        ls  $group/$sample/*_1.fastq* 2>/dev/null || ls $group/$sample/*_1.fq*  2>/dev/null || { echo Read file not found for $sample! Please make sure the fastq files are named as xxx_1.fastq, xxx_2.fastq or xxx_1.fq, xxx_2.fq;  exit 1; }
        read1=""
        read2=""
        
        for r1 in `ls $group/$sample/*_1.fastq* $group/$sample/*_1.fq* 2>/dev/null | xargs -n 1 basename`; do 
            #echo r1 is $r1
            readgroup=${r1%_*}
            echo working on readgroup: $readgroup
            r2=${r1%_*}_2${r1##*_1}
            #echo R2 is $r2
            mode=pe     
			[[ -f $group/$sample/$r2 ]] && r2="$group/$sample/$r2" || { echo -e "\n\n!!!Warning: read2 file '$r2' not exist, ignore this warning if you are working with single-end data\n\n"; r2=""; mode=any;  }
            out=skewer/$group.$sample.$readgroup
            r1=$group/$sample/$r1
           
            #@1,0,skewer,,sbatch -p short -t 2:0:0 -n 4 --mem 20G
            mkdir -p $out; skewer -t 4 -x $a -m $mode $r1 $r2 -o $out/SKEWER
            
            [ -z $r2 ] && read1="$read1,$pwd/$out/SKEWER-trimmed.fastq" || { read1="$read1,$pwd/$out/SKEWER-trimmed-pair1.fastq";  read2="$read2,$pwd/$out/SKEWER-trimmed-pair2.fastq"; }
            
            out=fastqc/$group.$sample.$readgroup; 
            
            #@2,0,fastqc,,sbatch -p short -t 2:0:0 -n 1 --mem 8G
            mkdir -p $out; fastqc -o $out $r1 $r2  
        done
        
        [ -z "$read2" ] && reads="-U ${read1#,}" || reads="-1 ${read1#,} -2 ${read2#,}"
        
        ############
        #echo reads: $reads  
        mkdir -p hisatOut/$group.$sample stringtieOut/$group.$sample ; cd hisatOut/$group.$sample; 
         
        #@3,1,hisatStringTie,,sbatch -n 4 -p short -t 12:0:0 --mem 40G
        rm aligns.sorted.bam.tmp* 2>/dev/null; hisat2 $index $reads --phred33 --mm -p 4 --dta  $splice | samtools view -Suh - |  samtools sort - -o aligns.sorted.bam && stringtie -G $gtf aligns.sorted.bam -o ../../stringtieOut/$group.$sample/stringtie.gtf
        
        #hisat2 $index $reads --phred33 --mm -p 4 --dta  $splice 2>$pwd/hisatOut/$group.$sample/hisat2.log | samtools view -Suh - |  samtools sort - -o aligns.sorted.bam; stringtie -G $gtf aligns.sorted.bam -o ../../stringtieOut/$group.$sample/stringtie.gtf
        
        gtfFiles="$gtfFiles stringtieOut/$group.$sample/stringtie.gtf"
        
       
        cd $pwd
        
        #exit  
    done 
done

#@4,3,stringtieMerge,,sbatch -n 4 -p short -t 2:0:0 --mem 40G
stringtie --merge -p 4 -o stringtieOut/merged.gtf -G $gtf $gtfFiles

for i in `ls -d stringtieOut/g*/`; do 
    echo workding on $i bam is: hisatOut/${i#stringtieOut/}aligns.sorted.bam 
    
    # see: https://ccb.jhu.edu/software/stringtie/index.shtml?t=manual
    
    #@5,4,stringtie2,  
    mkdir -p stringtieOut/sample${i#stringtieOut/}; stringtie -G stringtieOut/merged.gtf -b stringtieOut/sample${i#stringtieOut/}  -e -o $i/stringtie1.gtf hisatOut/${i#stringtieOut/}aligns.sorted.bam  
done

groups=${groups#,}
echo find final group membership for samples: $groups

#@6,5,ballgown,
Rscript ballgown.r $groups
