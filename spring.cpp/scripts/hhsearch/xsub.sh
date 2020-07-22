#!/bin/bash
TM=320:00

# wait
sleep 1

### job submission
cat<<JOB | qsub
#PBS -N spectrum_hhs
#PBS -q default
#PBS -d /tmp/
#PBS -o /tmp/inp.out
#PBS -e /tmp/inp.err
#PBS -l nodes=1:ppn=1,walltime=$TM
cd $1
$2 $3 $4 $5 $6 $7 $8 $9 $10
JOB
