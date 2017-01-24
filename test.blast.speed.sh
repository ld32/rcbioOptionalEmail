head -n 50000 all_faa.fa > in.fa
cp in.fa in1.fa

bsub -Is -n 4 -q interactive bash
on (clarinet002-062 )

module load seq/blast/ncbi-blast/2.2.26


(time blastp -db nr -query in.fa -num_threads 8 -out out.txt) > time 2>&1 

(time (mkdir -p /tmp/blast_db/ && mkdir -p /tmp/blast_in/ && mkdir -p /tmp/blast_out/ && rsync --include "nr*" --exclude "*" /hms/scratch1/shared_databases/blastdb/* /tmp/blast_db && rsync in.fa /tmp/blast_in && blastp -db /tmp/blast_db/nr -query /tmp/blast_in/in.fa -num_threads 8 -out /tmp/blast_out/out.txt && mv /tmp/blast_out/out.txt out1.txt)) > time1 2>&1 

module load seq/blast/ncbi-blast/2.2.26

(date && blastp -db nr -query in.fa -num_threads 4 -out out.txt && date) > run1.txt & 

(date && blast_wrapper.sh blastp -db nr -query in.fa -num_threads 4 -out out1.txt && date) > run2.txt &

 # local add a line
#add a line
# add another line
# add a line from remote master
