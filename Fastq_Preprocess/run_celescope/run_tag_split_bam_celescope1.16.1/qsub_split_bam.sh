#!/bin/bash
#-V:pass all environment variables to the job
#-N: job name
#-l: specify the amount of maximum memory required
#-d: working directory
work_dir=$PWD #default current working directory
echo $work_dir
#qsub -V -N qTest -l h_vmem=5G -d $work_dir test.sh
#qsub -V -N mt -l mem=10000MB -v "arg1=1000,arg2=444" -d $work_dir run1.sh
#qsub -V -N qs -l mem=10000MB,nodes=cu02:ppn=1 -d $work_dir run.sh
qsub -V -N qtag_split_bam -l mem=120000MB,nodes=cu02  -d $work_dir run_split_bam.sh
