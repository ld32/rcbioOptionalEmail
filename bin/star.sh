#!/bin/sh

usage() { echo -e "\nUsage : ${0##*/} <-l readLength, required> [-r speciesIndex (required if no -b. Such as: tair10,dm3,mm9,mm10,hg18 or hg19. Let us know if you need other references)] [-b starIndexWithPath(required if no -r, don't need this if -r is given)] [-f genomeFastaFile(required if no -r, don't need this if -r is given)] [-g gtfFileWithPath (optional, don't need this if -r is given)] "; exit 1;} 

while getopts ":r:b:f:g:l:" o; do
    case "${o}" in
        r)
            r=${OPTARG}
            ;;
        b)
            b=${OPTARG}
            ;;
        f)
            f=${OPTARG}
            ;;
        g)
            g=${OPTARG}
            ;;
        l)
            l=${OPTARG}
            ;;
       *)
            usage
            ;;
    esac
done

[ -z "$l" ] && echo Please give read length && usage 

module load  star/2.5.4a  samtools/1.3.1 picard/2.8.0

#SHARED_DATABASES=/n/shared_db
SHARED_DATABASES=/n/groups/shared_databases

# set the correct index and gene annotation file
if [ -z "${r}" ]; then
    if [ ! -z "$b" ]; then 
        [ -f $b/genomeParameters.txt ] || { echo -e "Error: \ngenome star index could not be found: $b"; usage; } 
        index=$b
    else
        usage
    fi
    [ -f ${f} ] || { echo -e "Error: \n1genome fasta file could not be found: ${f}"; usage;}  
    
    if [ ! -z "${g}" ]; then
       [ -f "${g}" ] || { echo -e "Error: \ngtf file could not be found: $g"; usage; } 
       gtf="-G $g"
    fi
    
else 
  case "$r" in
    "mm10")index="$SHARED_DATABASES/igenome/Mus_musculus/UCSC/mm10/Sequence/starIndex"
       gtf="-G $SHARED_DATABASES/igenome/Mus_musculus/UCSC/mm10/Annotation/Genes/genes.gtf"
    ;;

    "mm9")index="$SHARED_DATABASES/igenome/Mus_musculus/UCSC/mm9/Sequence/starIndex"
       gtf="-G $SHARED_DATABASES/igenome/Mus_musculus/UCSC/mm9/Annotation/Genes/genes.gtf" 
    ;;

    "hg18") index="$SHARED_DATABASES/igenome/Homo_sapiens/UCSC/hg18/Sequence/starIndex"
        gtf="-G $SHARED_DATABASES/igenome/Homo_sapiens/UCSC/hg18/Annotation/Genes/genes.gtf"
    ;;

    "hg19") index="$SHARED_DATABASES/igenome/Homo_sapiens/UCSC/hg19/Sequence/starIndex"
        gtf="-G $SHARED_DATABASES/igenome/Homo_sapiens/UCSC/hg19/Annotation/Genes/genes.gtf"
    ;;
  
     "dm3") index="$SHARED_DATABASES/igenome/Drosophila_melanogaster/UCSC/dm3/Sequence/starIndex"
       gtf="-G $SHARED_DATABASES/igenome/Drosophila_melanogaster/UCSC/dm3/Annotation/Genes/genes.gtf" 
    ;;
    
   # "dm3") index="$SHARED_DATABASES/dm3/star/2.5.4a"
   #        gtf="-G $index/genes.gtf"
   # ;;
    
    *)  echo "Index '$r' is not supported. Please email rchelp@hms.harvard.edu for help."; exit
    ;;
  esac
  fa=$index/genome.fa
  
  [ -f $fa ] || { echo genome fasta file not exist: $fa; usage; } 
fi



[ -d group2 ] || { echo group2 is not found. You need at least two groups to run this pipeline; exit 1; }



pwdhere=`pwd`

for group in `ls -v -d group*`; do
    echo working on group:  $group
   
    #loopStart,sample
    for sample in `ls -d $group/*  | xargs -n 1 basename`; do 
        echo working on sample: $sample
        inputsams=""
        ls  $group/$sample/*_1.fastq* >/dev/null 2>&1 || ls $group/$sample/*_1.fq*  >/dev/null 2>&1 || { echo Read file not found for $sample! Please make sure the fastq files are named as xxx_1.fastq, xxx_2.fastq or xxx_1.fq, xxx_2.fq;  exit 1; }
          
        #loopStart,readgroup
        for r1 in `ls $group/$sample/*_1.fastq* $group/$sample/*_1.fq* 2>/dev/null | xargs -n 1 basename`; do 
            readgroup=${r1##*/}
            readgroup=${readgroup%_*}
            r2=${r1%_*}_2${r1##*_1}
            echo working on readgroup: $readgroup
			r1=$pwdhere/$group/$sample/$r1 
            
            [[ -f $group/$sample/$r2 ]] && r2=$pwdhere/$group/$sample/$r2 || { echo -e "\n\n!!!Warning: read2 file '$r2' not exist, ignore this warning if you are working with single-end data\n\n"; r2=""; }
            
            zipcmd=""
            [[ $r1 == *gz ]] && zipcmd="--readFilesCommand zcat"  
            [[ $r1 == *bz2 ]] && zipcmd="--readFilesCommand bzcat"    
            
            mkdir -p starOut/$group$sample$readgroup-star.p1
            mkdir -p starOut/$group$sample$readgroup-star.ref
            mkdir -p starOut/$group$sample$readgroup-star.p2
            
            # first pass 
            #@1,0,star,index.gtf,sbatch -p short -n 4 -t 0-4:0 --mem 40G 
            cd $pwdhere/starOut/$group$sample$readgroup-star.p1 && STAR --genomeDir $index --readFilesIn $r1 $r2 --runThreadN 4 $zipcmd  && cd $pwdhere/starOut/$group$sample$readgroup-star.ref && STAR --runMode genomeGenerate --genomeDir ./ --genomeFastaFiles $fa --sjdbFileChrStartEnd $pwdhere/starOut/$group$sample$readgroup-star.p1/SJ.out.tab --sjdbOverhang $(($l - 1)) --runThreadN 4  &&  cd $pwdhere/starOut/$group$sample$readgroup-star.p2 && STAR --genomeDir $pwdhere/starOut/$group$sample$readgroup-star.ref --readFilesIn $r1 $r2 --runThreadN 4  $zipcmd --outSAMstrandField intronMotif    
            
            # # the awk command is to correct a issue: https://groups.google.com/forum/#!topic/rna-star/Ta1Z2u4bPfc
            # second pass
            ##@3,2,star3:star3 command 
            # && cd starOut/$group$sample$readgroup-star.p2 && STAR --genomeDir ../../starOut/$group$sample$readgroup-star.ref --readFilesIn $r1 $r2 --runThreadN $n  $zipcmd --outSAMstrandField intronMotif --genomeLoad LoadAndRemove && runAwk && cd $pwdhere          
                   
            # add readgroup
            ##@4,3,star4:star4 command 
            #java -Xmx4g -Djava.io.tmpdir=/tmp -jar $PICARD/picard.jar AddOrReplaceReadGroups I=$sample/star.p2.$readgroup/Aligned.out.sam O=$sample/star.p2.$readgroup/out.bam SO=coordinate  RGID=$sample.$readgroup RGLB=rglb  RGPU=rgpu RGPL=illumina RGSM=$sample
            
            inputsams="$inputsams INPUT=starOut/$group$sample$readgroup-star.p2/Aligned.out.sam"
            
        #loopEnd 	
        done
        
        
        cd $pwdhere           
        mkdir starOut/$group$sample
        
        #@2,1,merge,,sbatch -p short -n 1 -t 0-4:0 --mem 10G  
        java -Xmx12g -jar $PICARD/picard-2.8.0.jar MergeSamFiles VALIDATION_STRINGENCY=SILENT OUTPUT=starOut/$group$sample/accepted_hits.bam  $inputsams SORT_ORDER=coordinate && samtools index starOut/$group$sample/accepted_hits.bam
        
        #break
    done 
    #break
done


