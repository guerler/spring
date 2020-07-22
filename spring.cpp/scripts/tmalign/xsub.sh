#!/bin/bash
TM=1000:00

sleep 1

### job submission
cat<<JOB | qsub
#PBS -q default
#PBS -d /tmp/
#PBS -l nodes=1:ppn=1,walltime=$TM
#PBS -o /tmp/inp.out
#PBS -e /tmp/inp.err
#PBS -N spectrum_tmalign
cd $1
$2 $3 $4 $5 $6 $7 $8 $9 $10
JOB
