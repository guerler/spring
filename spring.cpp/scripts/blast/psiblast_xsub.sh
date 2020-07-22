#!/bin/bash
TM=60:00

sleep 1

### job submission  -num_iterations x
cat<<JOB | qsub
#PBS -q urgent
#PBS -d /tmp/
#PBS -N spectrum_blast
#PBS -l nodes=1:ppn=1,walltime=$TM
cd $1
./bin/psiblast -evalue 10 -num_alignments 100000 -query $2 -db $3 -out $4
JOB
