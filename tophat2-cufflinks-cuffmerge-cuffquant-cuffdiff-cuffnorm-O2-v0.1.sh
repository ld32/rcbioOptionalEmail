#!/bin/sh
echo Running $0 $@
xsub="sbatch -p long -t 24:0:0 -n 4"
stamp=$(date -d "today" +"%Y%m%d%H%M")
mkdir -p flag
if [ -f flag/alljobs.jid ]; then
    checkJobsSlurm  flag/alljobs.jid 
    [ $? == 1 ] && exit 0;
fi
cwd=`realpath ./flag`
[ -f flag/alljobs.jid ] && mv flag/alljobs.jid flag/alljobs.jid.old
printf "%-10s   %-20s   %-10s\n" job_id depend_on job_flag > flag/alljobs.jid
echo ---------------------------------------------------------
#!/bin/sh
usage() { echo -e "\nUsage : $0 [-q queueName(optional, default: long)] [-n numberOfCore(optional, default: 4)] [-W runTime(optional, default: 24:0:0)] [-r species_index (required if no -b. Such as: tair10,dm3,mm9,mm10,hg18 or hg19. Let us know if you need other references)] [-b bowtie2IndexWithPath(required if no -r, don't need this if -r is given)] [-f genomeFastaFile(required if no -r, don't need this if -r is given)] [-g gtfFileWithPath (optional, don't need this if -r is given)] [-t transcriptomeBotie2Index(optional, don't need this if -r is given)] [-p tophat parameters, default: --no-coverage-search ]"; exit 1;} 
while getopts ":q:n:W:r:b:f:g:t:p:" o; do
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
        b)
            b=${OPTARG}
            ;;
        f)
            f=${OPTARG}
            ;;
        g)
            g=${OPTARG}
            ;;
        t)
            t=${OPTARG}
            ;;
        p)  
            p=${OPTARG}
            ;;
        
        *)
            usage
            ;;
    esac
done
#echo 1 queue: $q core: $n run time: $W r: $r b: $b g: $g t: $t
[ -z "${q}" ] && q=long
[ -z "${n}" ] && n=4 
[ -z "${W}" ] && W=24:0:0
[ -z "${p}" ] && p="--no-coverage-search"
module load gcc/6.2.0 bowtie2/2.2.9 tophat/2.1.1 cufflinks/2.2.1 samtools/1.3.1
SHARED_DATABASES=/n/groups/shared_databases
xsub="sbatch -p $q -n $n -t $W"
# set the correct index and gene annotation file
if [ -z "${r}" ]; then
    if [ ! -z "$b" ]; then 
        bowtie2-inspect -n ${b} &> /dev/null || { echo -e "Error: \ngenome bowtie2 index could not be found: $b"; usage; } 
        index=$b
    else
        usage
    fi
    [ -f ${f} ] || { echo -e "Error: \ngenome fasta file could not be found: ${f}"; usage;}  
    
    if [ ! -z "${g}" ]; then
       [ -f "${g}" ] || { echo -e "Error: \ngtf file could not be found: $g"; usage; } 
       gtf="-G $g"
    fi
    if [ ! -z "$t" ]; then
        bowtie2-inspect -n ${t} &> /dev/null || { echo -e "Error: \ntranscriptome bowtie2 index could not be found: $t"; usage; } 
        trans="--transcriptome-index $t"   
    fi
else 
  case "$r" in
    "tair10")index="$SHARED_DATABASES/igenome/Arabidopsis_thaliana/NCBI/TAIR10/Sequence/Bowtie2Index/genome"
       gtf="-G $SHARED_DATABASES/igenome/Arabidopsis_thaliana/NCBI/TAIR10/Annotation/Genes/genes.gtf"
       trans="--transcriptome-index $SHARED_DATABASES/tophat2/Arabidopsis_thaliana/NCBI/TAIR10/Annotation/Genes/tophat2_trans/genes"       
    ;;
    "mm10")index="$SHARED_DATABASES/igenome/Mus_musculus/UCSC/mm10/Sequence/Bowtie2Index/genome"
       gtf="-G $SHARED_DATABASES/igenome/Mus_musculus/UCSC/mm10/Annotation/Genes/genes.gtf"
       trans="--transcriptome-index $SHARED_DATABASES/tophat2/Mus_musculus/UCSC/mm10/Annotation/Genes/tophat2_trans/genes"       
    ;;
    
    "mm9")index="$SHARED_DATABASES/igenome/Mus_musculus/UCSC/mm9/Sequence/Bowtie2Index/genome"
       gtf="-G $SHARED_DATABASES/igenome/Mus_musculus/UCSC/mm9/Annotation/Genes/genes.gtf" 
       trans="--transcriptome-index $SHARED_DATABASES/tophat2/Mus_musculus/UCSC/mm9/Annotation/Genes/tophat2_trans/genes" 
    ;;
    "dm3") index="$SHARED_DATABASES/igenome/Drosophila_melanogaster/UCSC/dm3/Sequence/Bowtie2Index/genome"
       gtf="-G $SHARED_DATABASES/igenome/Drosophila_melanogaster/UCSC/dm3/Annotation/Genes/genes.gtf" 
       trans="--transcriptome-index $SHARED_DATABASES/tophat2/Drosophila_melanogaster/UCSC/dm3/Annotation/Genes/tophat2_trans/genes" 
    ;;
    "dm3_2L_2R") index="$SHARED_DATABASES/igenome/Drosophila_melanogaster/UCSC/dm3_2L_2R/genome/bowtie2_index/genome"
       gtf="-G $SHARED_DATABASES/igenome/Drosophila_melanogaster/UCSC/dm3_2L_2R/Annotation/Genes/genes.gtf" 
       trans="--transcriptome-index $SHARED_DATABASES/tophat2/Drosophila_melanogaster/UCSC/dm3_2L_2R/Annotation/Genes/tophat2_trans/genes" 
    ;;	
    "dm5") index="$SHARED_DATABASES/igenome/Drosophila_melanogaster/NCBI/build5/Sequence/Bowtie2Index/genome"
       gtf="-G $SHARED_DATABASES/igenome/Drosophila_melanogaster/NCBI/build5/Annotation/Genes/genes.gtf"
       trans="--transcriptome-index $SHARED_DATABASES/igenome/Drosophila_melanogaster/NCBI/build5/trans/transcriptomeIndex"
    ;;
    "dm6") index="$SHARED_DATABASES/igenome/Drosophila_melanogaster/flybase/dmel_r6.05_FB2015_02/fasta/bowtie2Index/genome"
        gtf="-G $SHARED_DATABASES/igenome/Drosophila_melanogaster/flybase/dmel_r6.05_FB2015_02/gtf/genes.gtf"
        trans="--transcriptome-index $SHARED_DATABASES/igenome/Drosophila_melanogaster/flybase/dmel_r6.05_FB2015_02/transcriptomeIndex/genes"
     ;;    
    "hg18") index="$SHARED_DATABASES/igenome/Homo_sapiens/UCSC/hg18/Sequence/Bowtie2Index/genome"
        gtf="-G $SHARED_DATABASES/igenome/Homo_sapiens/UCSC/hg18/Annotation/Genes/genes.gtf"
        trans="--transcriptome-index $SHARED_DATABASES/tophat2/Homo_sapiens/UCSC/hg18/Annotation/Genes/tophat2_trans/genes"
    ;;
    "hg19") index="$SHARED_DATABASES/igenome/Homo_sapiens/UCSC/hg19/Sequence/Bowtie2Index/genome"
        gtf="-G $SHARED_DATABASES/igenome/Homo_sapiens/UCSC/hg19/Annotation/Genes/genes.gtf"
        trans="--transcriptome-index $SHARED_DATABASES/tophat2/Homo_sapiens/UCSC/hg19/Annotation/Genes/tophat2_trans/genes"
    ;;
     "GRCh37")index="$SHARED_DATABASES/igenome/Homo_sapiens/Ensembl/GRCh37/Sequence/Bowtie2Index/genome"
        gtf="-G $SHARED_DATABASES/igenome/Homo_sapiens/Ensembl/GRCh37/Annotation/Genes/genes.gtf"
        trans="--transcriptome-index $SHARED_DATABASES/tophat2/Homo_sapiens/Ensembl/GRCh37/Annotation/Genes/tophat2_trans/genes"
    ;;
    "GRCh38")index="$SHARED_DATABASES/igenome/Homo_sapiens/genecode/grch38.p7/genome/genome"
        gtf="-G $SHARED_DATABASES/igenome/Homo_sapiens/genecode/grch38.p7/annotation/genes.gtf"
        trans="--transcriptome-index $SHARED_DATABASES/igenome/Homo_sapiens/genecode/grch38.p7/transcriptIndex/transcriptome"
    ;;
    "chrpic1") index="$SHARED_DATABASES/igenome/painted_turtle/chrpic1/genome/bowtie2index/genome"
        gtf="-G $SHARED_DATABASES/igenome/painted_turtle/chrpic1/annotation/genes/genes.gtf"
        trans="--transcriptome-index $SHARED_DATABASES/tophat2/painted_turtle/chrpic1/annotation/genes/tophat2_trans/genes"
    ;;
	
    *)  echo "Index '$r' is not supported. Please email rchelp@hms.harvard.edu for help."; exit
    ;;
  esac
  fa=$index.fa
fi
[ -d group2 ] || { echo group2 is not found. You need at least two groups to run this pipeline; exit 1; }
rm manifest_file  -fr  2>/dev/null 
mkdir -p tophatOut
mkdir -p linkOut
mkdir -p quantOut
# go through sample groups
#loopStart,group
for group in `ls -d group*/|sed 's|[/]||g'`; do
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
        
        #echo tophat command: tophat -o tophatOut/$group$sample$readgroup -p $n $gtf  $trans $p $index ${read1#,} ${read2#,}
        
        tophats="$tophats $group$sample"
        
        #@1,0,tophat2:tophat2 command

        echo; echo step: 1, depends on: 0, job name: tophat2, flag: tophat2.$group.$sample  
        flag=1.0.tophat2.$group.$sample
        flag=${flag//\//_}
        deps=""
        [ -z "${deps// /}" ] && deps=null || deps=${deps// /.}
        id=$(sbatchRun $xsub -flag $deps $cwd $flag "rm -r tophatOut/$group$sample 2>/dev/null ; tophat -o tophatOut/$group$sample -p $n $gtf  $trans $p $index ${read1#,} ${read2#,} && samtools index tophatOut/$group$sample/accepted_hits.bam")
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
        
        #@2,1,cufflink: cufflink command

        echo; echo step: 2, depends on: 1, job name: cufflink, flag: cufflink.$group.$sample  
        flag=2.1.cufflink.$group.$sample
        flag=${flag//\//_}
        deps=""
        deps="${jobID[1]}"
        [ -z "${deps// /}" ] && deps=null || deps=${deps// /.}
        id=$(sbatchRun $xsub -flag $deps $cwd $flag "rm -r linkOut/$group$sample 2>/dev/null ; cufflinks -o linkOut/$group$sample -p $n -q tophatOut/$group$sample/accepted_hits.bam")
        if [ -z "$id" ]; then
            echo  job $flag is not submitted
            jobID[2]=""
        else
            alljobs="$alljobs $id"
            printf "%-10s  %-20s  %-10s\n" $id $deps $flag >> $cwd/alljobs.jid
            startNewLoop[1]="no"
            [ -z ${startNewLoop[2]} ] && jobIDs[2]="" && startNewLoop[2]="no" 
            jobID[2]=$id
            jobIDs[2]=${jobIDs[2]}.$id
        fi
        
        [ -z "$bams" ] && bams="quantOut/$group$sample/abundances.cxb" || bams="$bams,quantOut/$group$sample/abundances.cxb"
        
        [ -z "$samples" ] && samples=$sample || samples=$samples+$sample   
               
        echo linkOut/$group$sample/transcripts.gtf >> manifest_file
        
    #loopEnd 
    done 
    
    bams="$bams "
    samples="$samples,"
    
#loopEnd     
done
#@3,2,cuffmerge: cuffmerge command

echo; echo step: 3, depends on: 2, job name: cuffmerge, flag: cuffmerge  
flag=3.2.cuffmerge
flag=${flag//\//_}
deps=""
deps="${jobIDs[2]}"
[ -z "${deps// /}" ] && deps=null || deps=${deps// /.}
id=$(sbatchRun $xsub -flag $deps $cwd $flag "rm -r merge 2>/dev/null ; cuffmerge -o merge -p $n -g ${gtf#-G } -s $fa manifest_file")
if [ -z "$id" ]; then
    echo  job $flag is not submitted
    jobID[3]=""
else
    alljobs="$alljobs $id"
    printf "%-10s  %-20s  %-10s\n" $id $deps $flag >> $cwd/alljobs.jid
    startNewLoop[2]=""
    [ -z ${startNewLoop[3]} ] && jobIDs[3]="" && startNewLoop[3]="no" 
    jobID[3]=$id
    jobIDs[3]=${jobIDs[3]}.$id
fi
#loopStart,i
for i in $tophats; do 
    echo working $i;
   
    #@4,3,cuffquant

    echo; echo step: 4, depends on: 3, job name: cuffquant, flag: cuffquant.$i  
    flag=4.3.cuffquant.$i
    flag=${flag//\//_}
    deps=""
    deps="${jobIDs[3]}"
    [ -z "${deps// /}" ] && deps=null || deps=${deps// /.}
    id=$(sbatchRun $xsub -flag $deps $cwd $flag "rm quantOut/$i; cuffquant -o quantOut/$i merge/merged.gtf tophatOut/$i/accepted_hits.bam -p 4 -b $fa -u -q")
    if [ -z "$id" ]; then
        echo  job $flag is not submitted
        jobID[4]=""
    else
        alljobs="$alljobs $id"
        printf "%-10s  %-20s  %-10s\n" $id $deps $flag >> $cwd/alljobs.jid
        startNewLoop[3]=""
        [ -z ${startNewLoop[4]} ] && jobIDs[4]="" && startNewLoop[4]="no" 
        jobID[4]=$id
        jobIDs[4]=${jobIDs[4]}.$id
    fi
#loopEnd     
done    
bams=${bams// ,/ } ; samples=${samples//,+/,} ; samples=${samples%,}
#@5,4,cuffnorm: cuffnorm command

echo; echo step: 5, depends on: 4, job name: cuffnorm, flag: cuffnorm  
flag=5.4.cuffnorm
flag=${flag//\//_}
deps=""
deps="${jobIDs[4]}"
[ -z "${deps// /}" ] && deps=null || deps=${deps// /.}
id=$(sbatchRun $xsub -flag $deps $cwd $flag "rm norm -r 2>/dev/null ; cuffnorm -o norm merge/merged.gtf $bams -p $n -L $samples -q ")
if [ -z "$id" ]; then
    echo  job $flag is not submitted
    jobID[5]=""
else
    alljobs="$alljobs $id"
    printf "%-10s  %-20s  %-10s\n" $id $deps $flag >> $cwd/alljobs.jid
    startNewLoop[4]=""
    [ -z ${startNewLoop[5]} ] && jobIDs[5]="" && startNewLoop[5]="no" 
    jobID[5]=$id
    jobIDs[5]=${jobIDs[5]}.$id
fi
xsub="sbatch -p priority -n 8 -t 720:0:0"
#@6,4,cuffdiff: cuffdiff command

echo; echo step: 6, depends on: 4, job name: cuffdiff, flag: cuffdiff  
flag=6.4.cuffdiff
flag=${flag//\//_}
deps=""
deps="${jobIDs[4]}"
[ -z "${deps// /}" ] && deps=null || deps=${deps// /.}
id=$(sbatchRun $xsub -flag $deps $cwd $flag "rm diff -r 2>/dev/null ; cuffdiff -o diff merge/merged.gtf $bams -p $n -b $fa -L $samples -u -q ")
if [ -z "$id" ]; then
    echo  job $flag is not submitted
    jobID[6]=""
else
    alljobs="$alljobs $id"
    printf "%-10s  %-20s  %-10s\n" $id $deps $flag >> $cwd/alljobs.jid
    startNewLoop[4]=""
    [ -z ${startNewLoop[6]} ] && jobIDs[6]="" && startNewLoop[6]="no" 
    jobID[6]=$id
    jobIDs[6]=${jobIDs[6]}.$id
fi
cd $cwd/..
echo all submitted jobs: 
cat flag/alljobs.jid
echo ---------------------------------------------------------
[ -f flag/alljobs.jid.first ] || cp flag/alljobs.jid flag/alljobs.jid.first 
