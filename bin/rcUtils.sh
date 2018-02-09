#!/bin/sh

tophat2Usage() { 
                                                                          
    echo -e "\nUsage:  \n$(basename $1) \n[-r species_index (required if no -b. Such as: tair10,dm3,mm9,mm10,hg18,hg19,hg38,GRCh37 and GRCh38. Let us know if you need other references)] \n[-q queueName(optional, default: long)] \n[-n numberOfCore(optional, default: 4)] \n[-W runTime(optional, default: 24:0:0)] \n[-b bowtiew2IndexWithPath(required if no -r, don't need this if -r is given)] \n[-f genomeFastaFile(required if no -r, don't need this if -r is given)] \n[-g gtfFileWithPath (optional, don't need this if -r is given)] \n[-t transcriptomeBotie2Index(optional, don't need this if -r is given)] \n[-p tophat parameters, default: --no-coverage-search ] \nNote: we print each option in a new line to make it easir to read. We you run the software, please put all the options in the same line together with the command"; 
    exit 1;
} 

rsyncToTmp() {

  [ -z "$1" ] && { echo "Usage: $0 </tmp/indexPath> </tmp/gtfPath> </tmp/bowtieIndexPath> ... "; return; } 
  for v in $@; do
      echo Working to copy: $v, waiting lock...
      [[ $v == /tmp/* ]] || continue
      ls ${v#/tmp/rcbio}* >/dev/null 2>&1 || { echo Reference file or folder not exist: ${v#/tmp/rcbio}; continue; } 
      lockFile=/tmp/${v//\//-}
      
      while ! ( set -o noclobber; echo "$$" > "$lockFile") 2> /dev/null; do
        echo waiting for lock file: $lockFile
        sleep 30
      done
      echo Got lock: $lockFile. Copying data to: $v
      trap 'rm -f "$lockFile"; exit $?' INT TERM EXIT
      mkdir -p ${v%/*}
      rsync -aL ${v#/tmp/rcbio}* ${v%/*} 
      chmod -R a+rwx ${v%/*}
      find  ${v%/*} -type f -exec chmod -x '{}' \;
      echo Copying is done for $v
      rm -f "$lockFile"
      trap - INT TERM EXIT
  done
}
export -f rsyncToTmp

setPath() {
 
  [ -z "$1" ] && { echo "Usage: $0 index.gtf.bowtieIndex"; return; } 
  for v in ${1//./ }; do
      #echo Reset path for: $v
      vx=${!v}
      [[ $vx == /tmp/* ]] && continue
      #vx=`realpath $vx` 
      eval $v=/tmp/rcbio$vx
  done
  #echo new path: $gtf, $index
}
export -f setPath

setPathBack() {
 
  [ -z "$1" ] && { echo "Usage: $0 index.gtf.bowtieIndex"; return; } 
  for v in ${1//./ }; do
      #echo ResetBack path for: $v
      vx=${!v}
      [[ $vx != /tmp/* ]] && continue
      #vx=`realpath $vx` 
      eval $v=${vx#/tmp/rcbio}
  done
  #echo new path: $gtf, $index
}
export -f setPathBack

setTophat2Reference(){
  export SHARED_DATABASES=/n/groups/shared_databases
    case "$1" in
        "tair10")index="$SHARED_DATABASES/igenome/Arabidopsis_thaliana/NCBI/TAIR10/Sequence/Bowtie2Index/genome"
           gtf="$SHARED_DATABASES/igenome/Arabidopsis_thaliana/NCBI/TAIR10/Annotation/Genes/genes.gtf"
           trans="$SHARED_DATABASES/tophat2/Arabidopsis_thaliana/NCBI/TAIR10/Annotation/Genes/tophat2_trans/genes"       
        ;;
        "mm10")index="$SHARED_DATABASES/igenome/Mus_musculus/UCSC/mm10/Sequence/Bowtie2Index/genome"
           gtf="$SHARED_DATABASES/igenome/Mus_musculus/UCSC/mm10/Annotation/Genes/genes.gtf"
           trans="$SHARED_DATABASES/tophat2/Mus_musculus/UCSC/mm10/Annotation/Genes/tophat2_trans/genes"       
        ;;
        
        "mm9")index="$SHARED_DATABASES/igenome/Mus_musculus/UCSC/mm9/Sequence/Bowtie2Index/genome"
           gtf="$SHARED_DATABASES/igenome/Mus_musculus/UCSC/mm9/Annotation/Genes/genes.gtf" 
           trans="$SHARED_DATABASES/tophat2/Mus_musculus/UCSC/mm9/Annotation/Genes/tophat2_trans/genes" 
        ;;
        "dm3") index="$SHARED_DATABASES/igenome/Drosophila_melanogaster/UCSC/dm3/Sequence/Bowtie2Index/genome"
           gtf="$SHARED_DATABASES/igenome/Drosophila_melanogaster/UCSC/dm3/Annotation/Genes/genes.gtf" 
           trans="$SHARED_DATABASES/tophat2/Drosophila_melanogaster/UCSC/dm3/Annotation/Genes/tophat2_trans/genes" 
        ;;
        "dm3_2L_2R") index="$SHARED_DATABASES/igenome/Drosophila_melanogaster/UCSC/dm3_2L_2R/genome/bowtie2_index/genome"
           gtf="$SHARED_DATABASES/igenome/Drosophila_melanogaster/UCSC/dm3_2L_2R/Annotation/Genes/genes.gtf" 
           trans="$SHARED_DATABASES/tophat2/Drosophila_melanogaster/UCSC/dm3_2L_2R/Annotation/Genes/tophat2_trans/genes" 
        ;;	
        "dm5") index="$SHARED_DATABASES/igenome/Drosophila_melanogaster/NCBI/build5/Sequence/Bowtie2Index/genome"
           gtf="$SHARED_DATABASES/igenome/Drosophila_melanogaster/NCBI/build5/Annotation/Genes/genes.gtf"
           trans="$SHARED_DATABASES/igenome/Drosophila_melanogaster/NCBI/build5/trans/transcriptomeIndex"
        ;;
        "dm6") index="$SHARED_DATABASES/igenome/Drosophila_melanogaster/flybase/dmel_r6.05_FB2015_02/fasta/bowtie2Index/genome"
            gtf="$SHARED_DATABASES/igenome/Drosophila_melanogaster/flybase/dmel_r6.05_FB2015_02/gtf/genes.gtf"
            trans="$SHARED_DATABASES/igenome/Drosophila_melanogaster/flybase/dmel_r6.05_FB2015_02/transcriptomeIndex/genes"
         ;;    
        "hg18") index="$SHARED_DATABASES/igenome/Homo_sapiens/UCSC/hg18/Sequence/Bowtie2Index/genome"
            gtf="$SHARED_DATABASES/igenome/Homo_sapiens/UCSC/hg18/Annotation/Genes/genes.gtf"
            trans="$SHARED_DATABASES/tophat2/Homo_sapiens/UCSC/hg18/Annotation/Genes/tophat2_trans/genes"
        ;;
        "hg19") index="$SHARED_DATABASES/igenome/Homo_sapiens/UCSC/hg19/Sequence/Bowtie2Index/genome"
            gtf="$SHARED_DATABASES/igenome/Homo_sapiens/UCSC/hg19/Annotation/Genes/genes.gtf"
            trans="$SHARED_DATABASES/tophat2/Homo_sapiens/UCSC/hg19/Annotation/Genes/tophat2_trans/genes"
        ;;
         "hg38") index="$SHARED_DATABASES/igenome/Homo_sapiens/UCSC/hg38/Sequence/Bowtie2Index/genome"
            gtf="-G $SHARED_DATABASES/igenome/Homo_sapiens/UCSC/hg38/Annotation/Genes/genes.gtf"
            trans="--transcriptome-index $SHARED_DATABASES/tophat2/Homo_sapiens/UCSC/hg38/Annotation/Genes/tophat2_trans/genes"
        ;;
        "GRCh37")index="$SHARED_DATABASES/igenome/Homo_sapiens/Ensembl/GRCh37/Sequence/Bowtie2Index/genome"
            gtf="$SHARED_DATABASES/igenome/Homo_sapiens/Ensembl/GRCh37/Annotation/Genes/genes.gtf"
            trans="$SHARED_DATABASES/tophat2/Homo_sapiens/Ensembl/GRCh37/Annotation/Genes/tophat2_trans/genes"
        ;;

        "GRCh38")index="$SHARED_DATABASES/igenome/Homo_sapiens/genecode/grch38.p7/genome/genome"
            gtf="$SHARED_DATABASES/igenome/Homo_sapiens/genecode/grch38.p7/annotation/genes.gtf"
            trans="$SHARED_DATABASES/igenome/Homo_sapiens/genecode/grch38.p7/transcriptIndex/transcriptome"
        ;;

        "chrpic1") index="$SHARED_DATABASES/igenome/painted_turtle/chrpic1/genome/bowtie2index/genome"
            gtf="$SHARED_DATABASES/igenome/painted_turtle/chrpic1/annotation/genes/genes.gtf"
            trans="$SHARED_DATABASES/tophat2/painted_turtle/chrpic1/annotation/genes/tophat2_trans/genes"
        ;;
        
        *)  echo "Index '$r' is not supported. Please email rchelp@hms.harvard.edu for help."; 
        ;;
    esac

}
export -f setTophat2Reference


