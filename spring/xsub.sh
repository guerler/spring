#!/bin/bash
TM=120:00

sleep 1

### job submission
cat<<JOB | qsub
#PBS -q urgent
#PBS -d /tmp/
#PBS -l nodes=1:ppn=1,walltime=$TM
#PBS -o /tmp/inp.out
#PBS -e /tmp/inp.err
#PBS -N spring_rx
cd $1
$2 $3 $4 $5 $6 $7 $8 $9 $10
JOB
