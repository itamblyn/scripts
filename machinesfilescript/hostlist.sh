#!/bin/sh
#$ -S /bin/csh
#$ -cwd
#$ -N hostlist
#$ -j y
#$ -o hostlist.out
#$ -pe openmp 16 
#$ -l h_vmem=250M,h_rt=00:01:00
date
echo PE_HOSTFILE is $PE_HOSTFILE and contains:
echo -----------------------------
cat $PE_HOSTFILE
echo -----------------------------
echo HOSTNAME is $HOSTNAME
echo NHOSTS   is $NHOSTS
echo NSLOTS   is $NSLOTS
echo TMPDIR   is $TMPDIR
echo JOB_ID   is $JOB_ID
echo SGE_O_WORKDIR is $SGE_O_WORKDIR
# You might want to adjust this next path to your liking...
$SGE_O_WORKDIR/startmpi.sh $PE_HOSTFILE
echo $TMPDIR/machines contains:
echo -----------------------------
cat $TMPDIR/machines
echo -----------------------------
